#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run.py

双端测序数据处理主流程

流程概述：
1. 检查并生成必要文件（最长转录本、碱基转换文件、索引）
2. 三次mapping：STAR基因组正链、STAR基因组负链、bowtie2转录组
3. 合并三个BAM文件
4. 使用索引文件还原原始BAM并过滤
5. Pileup处理
6. 根据分子类型调用修饰位点
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional, Tuple

# 处理导入路径：支持直接运行和模块运行
# 检测是否直接运行（__name__ == "__main__" 或直接执行）
_script_dir = Path(__file__).parent
_scripts_dir = _script_dir / "scripts"

# 如果scripts目录不在路径中，添加它（用于直接运行）
if str(_scripts_dir) not in sys.path:
    sys.path.insert(0, str(_scripts_dir))
# 如果父目录不在路径中，添加它（用于模块导入）
if str(_script_dir.parent) not in sys.path:
    sys.path.insert(0, str(_script_dir.parent))

# 尝试绝对导入，失败则使用相对导入
try:
    from final_PE_pipeline.scripts.extract_longest_isoform import extract_longest_isoform
    from final_PE_pipeline.scripts.convert_fastq_fasta import change_fasta, change_fastq
    from final_PE_pipeline.scripts.build_index import build_index
    from final_PE_pipeline.scripts.mapping import mapping
    from final_PE_pipeline.scripts.extract_target_aln_dump_others import extract_target_aln_and_dump_others
    from final_PE_pipeline.scripts.transTogenome import change_trans_to_genome
    from final_PE_pipeline.scripts.merge_bam import merge_bams
    from final_PE_pipeline.scripts.return_origin_bam_and_filter import filter_and_return_origin_bam
    from final_PE_pipeline.scripts.pileup_PE import get_genome_wide_pileup_paired
    from final_PE_pipeline.scripts.call_site_FDR import call_site_FDR, get_modification_name
    from final_PE_pipeline.scripts.run_cmd_and_logging import setup_project_logger, get_logger
except ImportError:
    # 如果绝对导入失败，使用相对导入（直接运行模式）
    from scripts.extract_longest_isoform import extract_longest_isoform
    from scripts.convert_fastq_fasta import change_fasta, change_fastq
    from scripts.build_index import build_index
    from scripts.mapping import mapping
    from scripts.extract_target_aln_dump_others import extract_target_aln_and_dump_others
    from scripts.transTogenome import change_trans_to_genome
    from scripts.merge_bam import merge_bams
    from scripts.return_origin_bam_and_filter import filter_and_return_origin_bam
    from scripts.pileup_PE import get_genome_wide_pileup_paired
    from scripts.call_site_FDR import call_site_FDR, get_modification_name
    from scripts.run_cmd_and_logging import setup_project_logger, get_logger


def check_and_prepare_files_dna(
    genome_fa: Path,
    tooldir: Path,
    logger: logging.Logger,
    dryrun: bool = False
) -> Tuple[Path, Path, Path, Path]:
    """
    检查并生成DNA流程所需的文件（碱基转换文件、索引）
    
    返回:
        (genome_TC_GA_fa, genome_AG_CT_fa,
         genome_TC_GA_index, genome_AG_CT_index)
    """
    logger.info("=" * 60)
    logger.info("Step 1: Checking and preparing necessary files for DNA")
    logger.info("=" * 60)
    
    # 定义目录结构
    ref_dir = tooldir / "reference"
    index_dir = tooldir / "index"
    ref_dir.mkdir(parents=True, exist_ok=True)
    index_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. 检查并生成碱基转换文件
    # 基因组 TC_GA convert (T->C, G->A)
    genome_TC_GA_fa = ref_dir / f"{genome_fa.stem}_TC_GAconvert.fa"
    if not genome_TC_GA_fa.exists():
        if dryrun:
            logger.info(f"[DRYRUN] Would convert genome: TC_GA")
            logger.info(f"[DRYRUN]   Input: {genome_fa}")
            logger.info(f"[DRYRUN]   Output: {genome_TC_GA_fa}")
            logger.info(f"[DRYRUN]   Conversion rules: T:C, G:A")
        else:
            logger.info("Converting genome: TC_GA...")
            change_fasta(genome_fa, genome_TC_GA_fa, ['T:C', 'G:A'])
            logger.info(f"Converted genome (TC_GA): {genome_TC_GA_fa}")
    else:
        logger.info(f"Converted genome (TC_GA) already exists: {genome_TC_GA_fa}")
    
    # 基因组 AG_CT convert (A->G, C->T)
    genome_AG_CT_fa = ref_dir / f"{genome_fa.stem}_AG_CTconvert.fa"
    if not genome_AG_CT_fa.exists():
        if dryrun:
            logger.info(f"[DRYRUN] Would convert genome: AG_CT")
            logger.info(f"[DRYRUN]   Input: {genome_fa}")
            logger.info(f"[DRYRUN]   Output: {genome_AG_CT_fa}")
            logger.info(f"[DRYRUN]   Conversion rules: A:G, C:T")
        else:
            logger.info("Converting genome: AG_CT...")
            change_fasta(genome_fa, genome_AG_CT_fa, ['A:G', 'C:T'])
            logger.info(f"Converted genome (AG_CT): {genome_AG_CT_fa}")
    else:
        logger.info(f"Converted genome (AG_CT) already exists: {genome_AG_CT_fa}")
    
    # 2. 检查并生成索引
    # 基因组 TC_GA index (STAR)
    genome_TC_GA_index = index_dir / f"{genome_fa.stem}_TC_GA_STAR_index"
    if not genome_TC_GA_index.exists() or not (genome_TC_GA_index / "Genome").exists():
        if dryrun:
            logger.info(f"[DRYRUN] Would build STAR index for genome TC_GA")
            logger.info(f"[DRYRUN]   Genome: {genome_TC_GA_fa}")
            logger.info(f"[DRYRUN]   Output: {genome_TC_GA_index}")
            logger.info(f"[DRYRUN]   Threads: 16")
            logger.info(f"[DRYRUN]   Note: Building without GTF file")
        else:
            logger.info("Building STAR index for genome TC_GA (without GTF file)...")
            build_index(
                tool='STAR',
                genome_fasta=genome_TC_GA_fa,
                output_path=genome_TC_GA_index,
                threads=16,
                gtf_file=None,
                logger=logger
            )
            logger.info(f"STAR index (genome TC_GA): {genome_TC_GA_index}")
    else:
        logger.info(f"STAR index (genome TC_GA) already exists: {genome_TC_GA_index}")
    
    # 基因组 AG_CT index (STAR)
    genome_AG_CT_index = index_dir / f"{genome_fa.stem}_AG_CT_STAR_index"
    if not genome_AG_CT_index.exists() or not (genome_AG_CT_index / "Genome").exists():
        if dryrun:
            logger.info(f"[DRYRUN] Would build STAR index for genome AG_CT")
            logger.info(f"[DRYRUN]   Genome: {genome_AG_CT_fa}")
            logger.info(f"[DRYRUN]   Output: {genome_AG_CT_index}")
            logger.info(f"[DRYRUN]   Threads: 16")
            logger.info(f"[DRYRUN]   Note: Building without GTF file")
        else:
            logger.info("Building STAR index for genome AG_CT (without GTF file)...")
            build_index(
                tool='STAR',
                genome_fasta=genome_AG_CT_fa,
                output_path=genome_AG_CT_index,
                threads=16,
                gtf_file=None,
                logger=logger
            )
            logger.info(f"STAR index (genome AG_CT): {genome_AG_CT_index}")
    else:
        logger.info(f"STAR index (genome AG_CT) already exists: {genome_AG_CT_index}")
    
    return (genome_TC_GA_fa, genome_AG_CT_fa,
            genome_TC_GA_index, genome_AG_CT_index)


def check_and_prepare_files_rna(
    genome_fa: Path,
    trans_fa: Path,
    gtf_file: Path,
    tooldir: Path,
    logger: logging.Logger,
    dryrun: bool = False
) -> Tuple[Path, Path, Path, Path, Path, Path, Path]:
    """
    检查并生成RNA流程所需的文件（最长转录本、碱基转换文件、索引）
    
    返回:
        (genome_AG_CT_fa, genome_TC_GA_fa, longest_trans_AG_CT_fa,
         genome_AG_CT_index, genome_TC_GA_index, longest_trans_AG_CT_index,
         longest_trans_file)
    """
    logger.info("=" * 60)
    logger.info("Step 1: Checking and preparing necessary files for RNA")
    logger.info("=" * 60)
    
    # 定义目录结构
    ref_dir = tooldir / "reference"
    anno_dir = tooldir / "annotations"
    index_dir = tooldir / "index"
    ref_dir.mkdir(parents=True, exist_ok=True)
    anno_dir.mkdir(parents=True, exist_ok=True)
    index_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. 检查并生成最长转录本文件
    longest_trans_fa_raw = ref_dir / f"{genome_fa.stem}_longest_transcript.fa"
    longest_trans_file = anno_dir / f"{genome_fa.stem}_longest.trans"
    
    if not longest_trans_fa_raw.exists() or not longest_trans_file.exists():
        if dryrun:
            logger.info("[DRYRUN] Would extract longest isoforms")
            logger.info(f"[DRYRUN]   Input: {trans_fa}")
            logger.info(f"[DRYRUN]   GTF: {gtf_file}")
            logger.info(f"[DRYRUN]   Output FASTA: {longest_trans_fa_raw}")
            logger.info(f"[DRYRUN]   Output TSV: {longest_trans_file}")
        else:
            logger.info("Extracting longest isoforms...")
            extract_longest_isoform(
                input_fasta=trans_fa,
                gtf_path=gtf_file,
                output_fasta=longest_trans_fa_raw,
                output_trans_file=longest_trans_file
            )
            logger.info(f"Longest transcript FASTA: {longest_trans_fa_raw}")
            logger.info(f"Longest transcript info: {longest_trans_file}")
    else:
        logger.info(f"Longest transcript files already exist: {longest_trans_fa_raw}, {longest_trans_file}")
    
    # 2. 检查并生成碱基转换文件
    # 基因组 AG_CT convert
    genome_AG_CT_fa = ref_dir / f"{genome_fa.stem}_AG_CTconvert.fa"
    if not genome_AG_CT_fa.exists():
        if dryrun:
            logger.info(f"[DRYRUN] Would convert genome: AG_CT")
            logger.info(f"[DRYRUN]   Input: {genome_fa}")
            logger.info(f"[DRYRUN]   Output: {genome_AG_CT_fa}")
            logger.info(f"[DRYRUN]   Conversion rules: A:G, C:T")
        else:
            logger.info("Converting genome: AG_CT...")
            change_fasta(genome_fa, genome_AG_CT_fa, ['A:G', 'C:T'])
            logger.info(f"Converted genome (AG_CT): {genome_AG_CT_fa}")
    else:
        logger.info(f"Converted genome (AG_CT) already exists: {genome_AG_CT_fa}")
    
    # 基因组 TC_GA convert
    genome_TC_GA_fa = ref_dir / f"{genome_fa.stem}_TC_GAconvert.fa"
    if not genome_TC_GA_fa.exists():
        if dryrun:
            logger.info(f"[DRYRUN] Would convert genome: TC_GA")
            logger.info(f"[DRYRUN]   Input: {genome_fa}")
            logger.info(f"[DRYRUN]   Output: {genome_TC_GA_fa}")
            logger.info(f"[DRYRUN]   Conversion rules: T:C, G:A")
        else:
            logger.info("Converting genome: TC_GA...")
            change_fasta(genome_fa, genome_TC_GA_fa, ['T:C', 'G:A'])
            logger.info(f"Converted genome (TC_GA): {genome_TC_GA_fa}")
    else:
        logger.info(f"Converted genome (TC_GA) already exists: {genome_TC_GA_fa}")
    
    # 最长转录本 AG_CT convert
    longest_trans_AG_CT_fa = ref_dir / f"{genome_fa.stem}_longest_transcript_AG_CTconvert.fa"
    if not longest_trans_AG_CT_fa.exists():
        if dryrun:
            logger.info(f"[DRYRUN] Would convert longest transcript: AG_CT")
            logger.info(f"[DRYRUN]   Input: {longest_trans_fa_raw}")
            logger.info(f"[DRYRUN]   Output: {longest_trans_AG_CT_fa}")
            logger.info(f"[DRYRUN]   Conversion rules: A:G, C:T")
        else:
            logger.info("Converting longest transcript: AG_CT...")
            change_fasta(longest_trans_fa_raw, longest_trans_AG_CT_fa, ['A:G', 'C:T'])
            logger.info(f"Converted longest transcript (AG_CT): {longest_trans_AG_CT_fa}")
    else:
        logger.info(f"Converted longest transcript (AG_CT) already exists: {longest_trans_AG_CT_fa}")
    
    # 3. 检查并生成索引
    # 基因组 AG_CT index (STAR)
    genome_AG_CT_index = index_dir / f"{genome_fa.stem}_AG_CT_STAR_index"
    if not genome_AG_CT_index.exists() or not (genome_AG_CT_index / "Genome").exists():
        if dryrun:
            logger.info(f"[DRYRUN] Would build STAR index for genome AG_CT")
            logger.info(f"[DRYRUN]   Genome: {genome_AG_CT_fa}")
            logger.info(f"[DRYRUN]   Output: {genome_AG_CT_index}")
            logger.info(f"[DRYRUN]   Threads: 16")
            logger.info(f"[DRYRUN]   Note: Building without GTF file")
        else:
            logger.info("Building STAR index for genome AG_CT (without GTF file)...")
            build_index(
                tool='STAR',
                genome_fasta=genome_AG_CT_fa,
                output_path=genome_AG_CT_index,
                threads=16,
                gtf_file=None,  # 不使用GTF文件
                logger=logger
            )
            logger.info(f"STAR index (genome AG_CT): {genome_AG_CT_index}")
    else:
        logger.info(f"STAR index (genome AG_CT) already exists: {genome_AG_CT_index}")
    
    # 基因组 TC_GA index (STAR)
    genome_TC_GA_index = index_dir / f"{genome_fa.stem}_TC_GA_STAR_index"
    if not genome_TC_GA_index.exists() or not (genome_TC_GA_index / "Genome").exists():
        if dryrun:
            logger.info(f"[DRYRUN] Would build STAR index for genome TC_GA")
            logger.info(f"[DRYRUN]   Genome: {genome_TC_GA_fa}")
            logger.info(f"[DRYRUN]   Output: {genome_TC_GA_index}")
            logger.info(f"[DRYRUN]   Threads: 16")
            logger.info(f"[DRYRUN]   Note: Building without GTF file")
        else:
            logger.info("Building STAR index for genome TC_GA (without GTF file)...")
            build_index(
                tool='STAR',
                genome_fasta=genome_TC_GA_fa,
                output_path=genome_TC_GA_index,
                threads=16,
                gtf_file=None,  # 不使用GTF文件
                logger=logger
            )
            logger.info(f"STAR index (genome TC_GA): {genome_TC_GA_index}")
    else:
        logger.info(f"STAR index (genome TC_GA) already exists: {genome_TC_GA_index}")
    
    # 最长转录本 AG_CT index (bowtie2)
    longest_trans_AG_CT_index = index_dir / f"{genome_fa.stem}_longest_transcript_AG_CT_bowtie2_index"
    if not longest_trans_AG_CT_index.exists() or not (longest_trans_AG_CT_index.with_suffix('.1.bt2')).exists():
        if dryrun:
            logger.info(f"[DRYRUN] Would build bowtie2 index for longest transcript AG_CT")
            logger.info(f"[DRYRUN]   Transcriptome: {longest_trans_AG_CT_fa}")
            logger.info(f"[DRYRUN]   Output: {longest_trans_AG_CT_index}")
            logger.info(f"[DRYRUN]   Threads: 16")
        else:
            logger.info("Building bowtie2 index for longest transcript AG_CT...")
            build_index(
                tool='bowtie2',
                genome_fasta=longest_trans_AG_CT_fa,
                output_path=longest_trans_AG_CT_index,
                threads=16,
                logger=logger
            )
            logger.info(f"Bowtie2 index (longest transcript AG_CT): {longest_trans_AG_CT_index}")
    else:
        logger.info(f"Bowtie2 index (longest transcript AG_CT) already exists: {longest_trans_AG_CT_index}")
    
    return (genome_AG_CT_fa, genome_TC_GA_fa, longest_trans_AG_CT_fa,
            genome_AG_CT_index, genome_TC_GA_index, longest_trans_AG_CT_index,
            longest_trans_file)


def main_pipeline(
    r1_fastq: Path,
    r2_fastq: Path,
    genome_fa: Path,
    trans_fa: Path,
    gtf_file: Path,
    tooldir: Path,
    outputdir: Path,
    prefix: str,
    threads: int,
    type: str,
    dryrun: bool = False,
    logger: Optional[logging.Logger] = None,
    # Filter参数
    filter_target_bases: list = None,
    filter_cutoffs: list = None,
    # Call site参数
    min_pair_cov: int = 10,
    min_base_cov: int = 2,
    min_base_rate: float = 0.1,
    fdr_threshold: float = 0.05
):
    """
    主流程函数
    
    参数:
        r1_fastq: R1 FASTQ文件路径（原始，未转换）
        r2_fastq: R2 FASTQ文件路径（原始，未转换）
        genome_fa: 基因组FASTA文件路径
        trans_fa: 转录组FASTA文件路径
        gtf_file: GTF注释文件路径
        tooldir: 工具目录（项目根目录）
        outputdir: 输出目录
        prefix: 输出文件前缀
        threads: 线程数
        type: 分子类型，'DNA' 或 'RNA'
        dryrun: 是否仅打印不执行
        logger: logger对象（可选）
        filter_target_bases: 过滤目标碱基列表，默认['A', 'C']
        filter_cutoffs: 过滤阈值列表，默认[3, 3]
        min_pair_cov: 最小pair覆盖度，默认10
        min_base_cov: 最小碱基覆盖度，默认2
        min_base_rate: 最小碱基比率，默认0.1
        fdr_threshold: FDR阈值，默认0.05
    """
    if logger is None:
        logger = get_logger()
    
    type = type.upper()
    if type not in ['DNA', 'RNA']:
        raise ValueError(f"type must be 'DNA' or 'RNA', got '{type}'")
    
    # 根据type设置默认filter参数
    if filter_target_bases is None:
        if type == 'DNA':
            filter_target_bases = ['T', 'G']
        else:  # RNA
            filter_target_bases = ['A', 'C']
    if filter_cutoffs is None:
        # 所有cutoff使用相同的默认值
        filter_cutoffs = [3] * len(filter_target_bases)
    
    # 验证filter参数
    if len(filter_target_bases) != len(filter_cutoffs):
        raise ValueError(f"filter_target_bases和filter_cutoffs长度必须相同，"
                        f"got {len(filter_target_bases)} and {len(filter_cutoffs)}")
    
    # 记录所有参数到日志
    logger.info("=" * 60)
    logger.info("Pipeline Parameters")
    logger.info("=" * 60)
    logger.info(f"Input R1 FASTQ: {r1_fastq}")
    logger.info(f"Input R2 FASTQ: {r2_fastq}")
    logger.info(f"Genome FASTA: {genome_fa}")
    logger.info(f"Transcriptome FASTA: {trans_fa}")
    logger.info(f"GTF file: {gtf_file}")
    logger.info(f"Output directory: {outputdir}")
    logger.info(f"Prefix: {prefix}")
    logger.info(f"Threads: {threads}")
    logger.info(f"Type: {type}")
    logger.info(f"Filter target bases: {filter_target_bases}")
    logger.info(f"Filter cutoffs: {filter_cutoffs}")
    logger.info(f"Call site - min_pair_cov: {min_pair_cov}")
    logger.info(f"Call site - min_base_cov: {min_base_cov}")
    logger.info(f"Call site - min_base_rate: {min_base_rate}")
    logger.info(f"Call site - fdr_threshold: {fdr_threshold}")
    logger.info("=" * 60)
    
    # 设置输出目录
    outputdir.mkdir(parents=True, exist_ok=True)
    alignment_dir = outputdir / "alignment"
    alignment_dir.mkdir(parents=True, exist_ok=True)
    fastq_dir = outputdir / "fastq"
    fastq_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 0: 预处理FASTQ文件 - 转换并生成索引文件
    logger.info("=" * 60)
    logger.info("Step 0: Preprocessing FASTQ files - Converting and generating index files")
    logger.info("=" * 60)
    
    # R1: A:G, C:T 转换
    r1_converted_fq = fastq_dir / f"{prefix}_r1_converted.fq.gz"
    r1_index_file = fastq_dir / f"{prefix}_r1_index_A-G_C-T.txt"
    
    if dryrun:
        logger.info(f"[DRYRUN] Would convert R1 FASTQ and generate index")
        logger.info(f"[DRYRUN]   Input: {r1_fastq}")
        logger.info(f"[DRYRUN]   Conversion rules: A:G, C:T")
        logger.info(f"[DRYRUN]   Output converted FASTQ: {r1_converted_fq}")
        logger.info(f"[DRYRUN]   Output index: {r1_index_file}")
    else:
        logger.info("Converting R1 FASTQ: A:G, C:T...")
        r1_converted, r1_idx = change_fastq(
            input_fastq=r1_fastq,
            output_dir=fastq_dir,
            convert_list=['A:G', 'C:T']
        )
        # 确保使用正确的文件名
        if r1_converted != r1_converted_fq:
            if r1_converted.exists():
                r1_converted.rename(r1_converted_fq)
        if r1_idx != r1_index_file:
            if r1_idx.exists():
                r1_idx.rename(r1_index_file)
            else:
                logger.warning(f"Expected index file {r1_index_file}, but got {r1_idx}")
                r1_index_file = r1_idx
        logger.info(f"R1 converted FASTQ: {r1_converted_fq}")
        logger.info(f"R1 index file: {r1_index_file}")
    
    # R2: T:C, G:A 转换
    r2_converted_fq = fastq_dir / f"{prefix}_r2_converted.fq.gz"
    r2_index_file = fastq_dir / f"{prefix}_r2_index_T-C_G-A.txt"
    
    if dryrun:
        logger.info(f"[DRYRUN] Would convert R2 FASTQ and generate index")
        logger.info(f"[DRYRUN]   Input: {r2_fastq}")
        logger.info(f"[DRYRUN]   Conversion rules: T:C, G:A")
        logger.info(f"[DRYRUN]   Output converted FASTQ: {r2_converted_fq}")
        logger.info(f"[DRYRUN]   Output index: {r2_index_file}")
    else:
        logger.info("Converting R2 FASTQ: T:C, G:A...")
        r2_converted, r2_idx = change_fastq(
            input_fastq=r2_fastq,
            output_dir=fastq_dir,
            convert_list=['T:C', 'G:A']
        )
        # 确保使用正确的文件名
        if r2_converted != r2_converted_fq:
            if r2_converted.exists():
                r2_converted.rename(r2_converted_fq)
        if r2_idx != r2_index_file:
            if r2_idx.exists():
                r2_idx.rename(r2_index_file)
            else:
                logger.warning(f"Expected index file {r2_index_file}, but got {r2_idx}")
                r2_index_file = r2_idx
        logger.info(f"R2 converted FASTQ: {r2_converted_fq}")
        logger.info(f"R2 index file: {r2_index_file}")
    
    # 使用转换后的FASTQ进行后续流程
    r1_fq_for_mapping = r1_converted_fq
    r2_fq_for_mapping = r2_converted_fq
    
    # Step 1: 检查并准备必要文件（根据type）
    if dryrun:
        logger.info("=" * 60)
        logger.info("DRYRUN MODE: No actual operations will be performed")
        logger.info("=" * 60)
    
    if type == 'DNA':
        # DNA流程：只需要基因组转换文件和索引
        (genome_TC_GA_fa, genome_AG_CT_fa,
         genome_TC_GA_index, genome_AG_CT_index) = check_and_prepare_files_dna(
            genome_fa, tooldir, logger, dryrun
        )
        longest_trans_file = None
        longest_trans_AG_CT_fa = None
        longest_trans_AG_CT_index = None
    else:  # RNA
        # RNA流程：需要基因组、转录组转换文件和索引
        (genome_AG_CT_fa, genome_TC_GA_fa, longest_trans_AG_CT_fa,
         genome_AG_CT_index, genome_TC_GA_index, longest_trans_AG_CT_index,
         longest_trans_file) = check_and_prepare_files_rna(
            genome_fa, trans_fa, gtf_file, tooldir, logger, dryrun
        )
    
    # Step 2: 开始mapping流程（根据type）
    if type == 'DNA':
        # DNA流程：两步mapping到基因组
        # Step 2a: 第一次mapping - STAR基因组TC_GA正链
        logger.info("=" * 60)
        logger.info("Step 2a: First mapping - STAR genome TC_GA (R1 positive strand) - DNA")
        logger.info("=" * 60)
        
        bam1_prefix = alignment_dir / f"{prefix}.star_genome_TC_GA"
        bam1_raw = alignment_dir / f"{prefix}.star_genome_TC_GAAligned.sortedByCoord.out.bam"
        bam1_filtered = alignment_dir / f"{prefix}.star_genome_TC_GA.filtered.bam"
        remaining_r1_1 = alignment_dir / f"{prefix}.remaining_1.r1.fq.gz"
        remaining_r2_1 = alignment_dir / f"{prefix}.remaining_1.r2.fq.gz"
        
        if dryrun:
            logger.info(f"[DRYRUN] Would run STAR mapping")
            logger.info(f"[DRYRUN]   R1: {r1_fq_for_mapping}")
            logger.info(f"[DRYRUN]   R2: {r2_fq_for_mapping}")
            logger.info(f"[DRYRUN]   Index: {genome_TC_GA_index}")
            logger.info(f"[DRYRUN]   Output: {bam1_raw}")
            logger.info(f"[DRYRUN]   Type: {type}, Threads: {threads}")
            logger.info(f"[DRYRUN] Would extract R1 positive strand pairs")
            logger.info(f"[DRYRUN]   Output BAM: {bam1_filtered}")
            logger.info(f"[DRYRUN]   Remaining R1: {remaining_r1_1}")
            logger.info(f"[DRYRUN]   Remaining R2: {remaining_r2_1}")
            repaired_r1_1 = remaining_r1_1
            repaired_r2_1 = remaining_r2_1
        else:
            bam1_raw = mapping(
                tool='STAR',
                r1_fastq=r1_fq_for_mapping,
                r2_fastq=r2_fq_for_mapping,
                index_path=genome_TC_GA_index,
                output_path=bam1_prefix,
                type=type,
                threads=threads,
                logger=logger
            )
            
            _, repaired_r1_1, repaired_r2_1 = extract_target_aln_and_dump_others(
                bam_path=bam1_raw,
                output_bam=bam1_filtered,
                output_r1_fastq=remaining_r1_1,
                output_r2_fastq=remaining_r2_1,
                r1_strand='plus',
                logger=logger
            )
        
        # Step 2b: 第二次mapping - STAR基因组AG_CT负链
        logger.info("=" * 60)
        logger.info("Step 2b: Second mapping - STAR genome AG_CT (R1 negative strand) - DNA")
        logger.info("=" * 60)
        
        bam2_prefix = alignment_dir / f"{prefix}.star_genome_AG_CT"
        bam2_raw = alignment_dir / f"{prefix}.star_genome_AG_CTAligned.sortedByCoord.out.bam"
        bam2_filtered = alignment_dir / f"{prefix}.star_genome_AG_CT.filtered.bam"
        
        if dryrun:
            logger.info(f"[DRYRUN] Would run STAR mapping")
            logger.info(f"[DRYRUN]   R1: {repaired_r1_1}")
            logger.info(f"[DRYRUN]   R2: {repaired_r2_1}")
            logger.info(f"[DRYRUN]   Index: {genome_AG_CT_index}")
            logger.info(f"[DRYRUN]   Output: {bam2_raw}")
            logger.info(f"[DRYRUN]   Type: {type}, Threads: {threads}")
            logger.info(f"[DRYRUN] Would extract R1 negative strand pairs")
            logger.info(f"[DRYRUN]   Output BAM: {bam2_filtered}")
        else:
            bam2_raw = mapping(
                tool='STAR',
                r1_fastq=repaired_r1_1,
                r2_fastq=repaired_r2_1,
                index_path=genome_AG_CT_index,
                output_path=bam2_prefix,
                type=type,
                threads=threads,
                logger=logger
            )
            
            _, _, _ = extract_target_aln_and_dump_others(
                bam_path=bam2_raw,
                output_bam=bam2_filtered,
                output_r1_fastq=Path("/dev/null"),  # DNA流程不需要剩余reads
                output_r2_fastq=Path("/dev/null"),
                r1_strand='minus',
                logger=logger
            )
        
        # Step 2c: 合并两个BAM文件（DNA流程）
        logger.info("=" * 60)
        logger.info("Step 2c: Merging two BAM files - DNA")
        logger.info("=" * 60)
        
        merged_bam = alignment_dir / f"{prefix}.merged.sorted.bam"
        if dryrun:
            logger.info(f"[DRYRUN] Would merge two BAM files")
            logger.info(f"[DRYRUN]   BAM 1: {bam1_filtered}")
            logger.info(f"[DRYRUN]   BAM 2: {bam2_filtered}")
            logger.info(f"[DRYRUN]   Output: {merged_bam}")
        else:
            merge_bams(
                bams=[bam1_filtered, bam2_filtered],
                outputdir=alignment_dir,
                prefix=prefix,
                logger=logger
            )
    
    else:  # RNA流程
        # RNA流程：三步mapping（两次基因组+一次转录组）
        # Step 2a: 第一次mapping - STAR基因组AG_CT正链
        logger.info("=" * 60)
        logger.info("Step 2a: First mapping - STAR genome AG_CT (R1 positive strand) - RNA")
        logger.info("=" * 60)
        
        bam1_prefix = alignment_dir / f"{prefix}.star_genome_AG_CT"
        bam1_raw = alignment_dir / f"{prefix}.star_genome_AG_CTAligned.sortedByCoord.out.bam"
        bam1_filtered = alignment_dir / f"{prefix}.star_genome_AG_CT.filtered.bam"
        remaining_r1_1 = alignment_dir / f"{prefix}.remaining_1.r1.fq.gz"
        remaining_r2_1 = alignment_dir / f"{prefix}.remaining_1.r2.fq.gz"
        
        if dryrun:
            logger.info(f"[DRYRUN] Would run STAR mapping")
            logger.info(f"[DRYRUN]   R1: {r1_fq_for_mapping}")
            logger.info(f"[DRYRUN]   R2: {r2_fq_for_mapping}")
            logger.info(f"[DRYRUN]   Index: {genome_AG_CT_index}")
            logger.info(f"[DRYRUN]   Output: {bam1_raw}")
            logger.info(f"[DRYRUN]   Type: {type}, Threads: {threads}")
            logger.info(f"[DRYRUN] Would extract R1 positive strand pairs")
            logger.info(f"[DRYRUN]   Output BAM: {bam1_filtered}")
            logger.info(f"[DRYRUN]   Remaining R1: {remaining_r1_1}")
            logger.info(f"[DRYRUN]   Remaining R2: {remaining_r2_1}")
            repaired_r1_1 = remaining_r1_1
            repaired_r2_1 = remaining_r2_1
        else:
            bam1_raw = mapping(
                tool='STAR',
                r1_fastq=r1_fq_for_mapping,
                r2_fastq=r2_fq_for_mapping,
                index_path=genome_AG_CT_index,
                output_path=bam1_prefix,
                type=type,
                threads=threads,
                logger=logger
            )
            
            _, repaired_r1_1, repaired_r2_1 = extract_target_aln_and_dump_others(
                bam_path=bam1_raw,
                output_bam=bam1_filtered,
                output_r1_fastq=remaining_r1_1,
                output_r2_fastq=remaining_r2_1,
                r1_strand='plus',
                logger=logger
            )
        
        # Step 2b: 第二次mapping - STAR基因组TC_GA负链
        logger.info("=" * 60)
        logger.info("Step 2b: Second mapping - STAR genome TC_GA (R1 negative strand) - RNA")
        logger.info("=" * 60)
        
        bam2_prefix = alignment_dir / f"{prefix}.star_genome_TC_GA"
        bam2_raw = alignment_dir / f"{prefix}.star_genome_TC_GAAligned.sortedByCoord.out.bam"
        bam2_filtered = alignment_dir / f"{prefix}.star_genome_TC_GA.filtered.bam"
        remaining_r1_2 = alignment_dir / f"{prefix}.remaining_2.r1.fq.gz"
        remaining_r2_2 = alignment_dir / f"{prefix}.remaining_2.r2.fq.gz"
        
        if dryrun:
            logger.info(f"[DRYRUN] Would run STAR mapping")
            logger.info(f"[DRYRUN]   R1: {repaired_r1_1}")
            logger.info(f"[DRYRUN]   R2: {repaired_r2_1}")
            logger.info(f"[DRYRUN]   Index: {genome_TC_GA_index}")
            logger.info(f"[DRYRUN]   Output: {bam2_raw}")
            logger.info(f"[DRYRUN]   Type: {type}, Threads: {threads}")
            logger.info(f"[DRYRUN] Would extract R1 negative strand pairs")
            logger.info(f"[DRYRUN]   Output BAM: {bam2_filtered}")
            logger.info(f"[DRYRUN]   Remaining R1: {remaining_r1_2}")
            logger.info(f"[DRYRUN]   Remaining R2: {remaining_r2_2}")
            repaired_r1_2 = remaining_r1_2
            repaired_r2_2 = remaining_r2_2
        else:
            bam2_raw = mapping(
                tool='STAR',
                r1_fastq=repaired_r1_1,
                r2_fastq=repaired_r2_1,
                index_path=genome_TC_GA_index,
                output_path=bam2_prefix,
                type=type,
                threads=threads,
                logger=logger
            )
            
            _, repaired_r1_2, repaired_r2_2 = extract_target_aln_and_dump_others(
                bam_path=bam2_raw,
                output_bam=bam2_filtered,
                output_r1_fastq=remaining_r1_2,
                output_r2_fastq=remaining_r2_2,
                r1_strand='minus',
                logger=logger
            )
        
        # Step 2c: 第三次mapping - bowtie2转录组AG_CT正链
        logger.info("=" * 60)
        logger.info("Step 2c: Third mapping - bowtie2 transcriptome AG_CT (R1 positive strand) - RNA")
        logger.info("=" * 60)
        
        bam3_sam = alignment_dir / f"{prefix}.bowtie2_transcriptome_AG_CT.sam"
        bam3_raw = alignment_dir / f"{prefix}.bowtie2_transcriptome_AG_CT.bam"
        bam3_filtered = alignment_dir / f"{prefix}.bowtie2_transcriptome_AG_CT.filtered.bam"
        remaining_r1_3 = alignment_dir / f"{prefix}.remaining_3.r1.fq.gz"
        remaining_r2_3 = alignment_dir / f"{prefix}.remaining_3.r2.fq.gz"
        
        if dryrun:
            logger.info(f"[DRYRUN] Would run bowtie2 mapping")
            logger.info(f"[DRYRUN]   R1: {repaired_r1_2}")
            logger.info(f"[DRYRUN]   R2: {repaired_r2_2}")
            logger.info(f"[DRYRUN]   Index: {longest_trans_AG_CT_index}")
            logger.info(f"[DRYRUN]   Output: {bam3_raw}")
            logger.info(f"[DRYRUN]   Type: {type}, Threads: {threads}")
            logger.info(f"[DRYRUN] Would extract R1 positive strand pairs")
            logger.info(f"[DRYRUN]   Output BAM: {bam3_filtered}")
            logger.info(f"[DRYRUN]   Remaining R1: {remaining_r1_3}")
            logger.info(f"[DRYRUN]   Remaining R2: {remaining_r2_3}")
        else:
            bam3_raw = mapping(
                tool='bowtie2',
                r1_fastq=repaired_r1_2,
                r2_fastq=repaired_r2_2,
                index_path=longest_trans_AG_CT_index,
                output_path=bam3_sam,
                type=type,
                threads=threads,
                logger=logger
            )
            
            _, _, _ = extract_target_aln_and_dump_others(
                bam_path=bam3_raw,
                output_bam=bam3_filtered,
                output_r1_fastq=remaining_r1_3,
                output_r2_fastq=remaining_r2_3,
                r1_strand='plus',
                logger=logger
            )
        
        # Step 2d: 转录组到基因组坐标转换
        logger.info("=" * 60)
        logger.info("Step 2d: Converting transcriptome BAM to genome coordinates - RNA")
        logger.info("=" * 60)
        
        bam3_genome = alignment_dir / f"{prefix}.bowtie2_transcriptome_AG_CT.genome.bam"
        if dryrun:
            logger.info(f"[DRYRUN] Would convert transcriptome BAM to genome coordinates")
            logger.info(f"[DRYRUN]   Input BAM: {bam3_filtered}")
            logger.info(f"[DRYRUN]   Genome FASTA: {genome_AG_CT_fa}")
            logger.info(f"[DRYRUN]   Transcript info: {longest_trans_file}")
            logger.info(f"[DRYRUN]   Output BAM: {bam3_genome}")
        else:
            change_trans_to_genome(
                trans_bam=bam3_filtered,
                target_genome_fa=genome_AG_CT_fa,
                longest_transfile=longest_trans_file,
                outputbam=bam3_genome
            )
        
        # Step 2e: 合并三个BAM文件（RNA流程）
        logger.info("=" * 60)
        logger.info("Step 2e: Merging three BAM files - RNA")
        logger.info("=" * 60)
        
        merged_bam = alignment_dir / f"{prefix}.merged.sorted.bam"
        if dryrun:
            logger.info(f"[DRYRUN] Would merge three BAM files")
            logger.info(f"[DRYRUN]   BAM 1: {bam1_filtered}")
            logger.info(f"[DRYRUN]   BAM 2: {bam2_filtered}")
            logger.info(f"[DRYRUN]   BAM 3: {bam3_genome}")
            logger.info(f"[DRYRUN]   Output: {merged_bam}")
        else:
            merge_bams(
                bams=[bam1_filtered, bam2_filtered, bam3_genome],
                outputdir=alignment_dir,
                prefix=prefix,
                logger=logger
            )
    
    # Step 7: 还原原始BAM并过滤
    logger.info("=" * 60)
    logger.info("Step 8: Restoring original BAM and filtering")
    logger.info("=" * 60)
    
    origin_bam = alignment_dir / f"{prefix}.origin.bam"
    filter_bam = alignment_dir / f"{prefix}.filtered.bam"
    
    # 对每个target_base进行过滤
    # 首先恢复原始BAM（只在第一次）
    current_bam = merged_bam
    is_first_filter = True
    
    for i, (target_base, cutoff) in enumerate(zip(filter_target_bases, filter_cutoffs)):
        step_filter_bam = alignment_dir / f"{prefix}.filtered_{target_base}_{cutoff}.bam"
        
        if dryrun:
            logger.info(f"[DRYRUN] Would filter BAM with target_base={target_base}, cutoff={cutoff}")
            logger.info(f"[DRYRUN]   Input BAM: {current_bam}")
            logger.info(f"[DRYRUN]   Output BAM: {step_filter_bam}")
            if is_first_filter:
                logger.info(f"[DRYRUN]   Origin BAM: {origin_bam} (will be generated)")
        else:
            logger.info(f"Filtering BAM (step {i+1}/{len(filter_target_bases)}): target_base={target_base}, cutoff={cutoff}")
            # 第一次过滤时生成origin_bam，后续只过滤
            if is_first_filter:
                origin_bam_result, step_filter_bam_result = filter_and_return_origin_bam(
                    bam_path=current_bam,
                    r1_index_file=r1_index_file,
                    r2_index_file=r2_index_file,
                    origin_bam=origin_bam,
                    filter_bam=step_filter_bam,
                    target_base=target_base,
                    cutoff=cutoff,
                    logger=logger
                )
                is_first_filter = False
            else:
                # 后续过滤：只过滤，不恢复（因为已经是原始BAM了）
                # 但我们需要重新恢复，因为filter_and_return_origin_bam会同时恢复和过滤
                # 所以每次都需要恢复
                _, step_filter_bam_result = filter_and_return_origin_bam(
                    bam_path=current_bam,
                    r1_index_file=r1_index_file,
                    r2_index_file=r2_index_file,
                    origin_bam=current_bam,  # 使用当前BAM作为origin（已经是恢复后的）
                    filter_bam=step_filter_bam,
                    target_base=target_base,
                    cutoff=cutoff,
                    logger=logger
                )
            current_bam = step_filter_bam_result
    
    final_bam = current_bam
    if dryrun:
        logger.info(f"[DRYRUN] Final filtered BAM: {final_bam}")
    
    # Step 9: Pileup处理
    logger.info("=" * 60)
    logger.info("Step 9: Generating pileup")
    logger.info("=" * 60)
    
    pileup_file = outputdir / f"{prefix}.pileup"
    if dryrun:
        logger.info(f"[DRYRUN] Would generate pileup")
        logger.info(f"[DRYRUN]   Input BAM: {final_bam}")
        logger.info(f"[DRYRUN]   Reference FASTA: {genome_fa} (原始基因组，未转换)")
        logger.info(f"[DRYRUN]   Output: {pileup_file}")
        logger.info(f"[DRYRUN]   Threads: {threads}")
        logger.info(f"[DRYRUN]   Allowed ref bases: A, C")
    else:
        # 使用原始基因组FASTA，而不是转换后的
        # 根据type确定主要read和允许的ref bases
        if type == 'DNA':
            primary_read = 'R2'  # DNA按R2计数
            allowed_ref_bases = {'A', 'C'}  # DNA统计A和C
        else:  # RNA
            primary_read = 'R1'  # RNA按R1计数
            allowed_ref_bases = {'A', 'C'}  # RNA统计A和C
        
        get_genome_wide_pileup_paired(
            bam_path=final_bam,
            fasta_path=genome_fa,  # 使用原始基因组
            outputdir=outputdir,
            thread=threads,
            prefix=prefix,
            tmp_dir=None,
            chunk_size=500000,
            max_depth=100000,
            allowed_ref_bases=allowed_ref_bases,
            primary_read=primary_read
        )
    
    # Step 10: 调用修饰位点
    logger.info("=" * 60)
    logger.info("Step 10: Calling modification sites")
    logger.info("=" * 60)
    
    # 根据type调用不同的修饰位点
    # DNA: 6mA和5mC
    # RNA: m6A和m5C
    
    # 调用A位点
    mod_name_A = get_modification_name('A', type)
    if dryrun:
        logger.info(f"[DRYRUN] Would call A modification sites")
        logger.info(f"[DRYRUN]   Pileup: {pileup_file}")
        logger.info(f"[DRYRUN]   Molecule type: {type}")
        logger.info(f"[DRYRUN]   Modification: {mod_name_A}")
        logger.info(f"[DRYRUN]   min_pair_cov: {min_pair_cov}")
        logger.info(f"[DRYRUN]   min_base_cov: {min_base_cov}")
        logger.info(f"[DRYRUN]   min_base_rate: {min_base_rate}")
        logger.info(f"[DRYRUN]   fdr_threshold: {fdr_threshold}")
    else:
        call_site_FDR(
            pileup_file=pileup_file,
            output_dir=outputdir,
            prefix=prefix,
            base_type='A',
            molecule_type=type,
            min_pair_cov=min_pair_cov,
            min_base_cov=min_base_cov,
            min_base_rate=min_base_rate,
            fdr_threshold=fdr_threshold,
            logger=logger
        )
    
    # 调用C位点
    mod_name_C = get_modification_name('C', type)
    if dryrun:
        logger.info(f"[DRYRUN] Would call C modification sites")
        logger.info(f"[DRYRUN]   Pileup: {pileup_file}")
        logger.info(f"[DRYRUN]   Molecule type: {type}")
        logger.info(f"[DRYRUN]   Modification: {mod_name_C}")
        logger.info(f"[DRYRUN]   min_pair_cov: {min_pair_cov}")
        logger.info(f"[DRYRUN]   min_base_cov: {min_base_cov}")
        logger.info(f"[DRYRUN]   min_base_rate: {min_base_rate}")
        logger.info(f"[DRYRUN]   fdr_threshold: {fdr_threshold}")
    else:
        call_site_FDR(
            pileup_file=pileup_file,
            output_dir=outputdir,
            prefix=prefix,
            base_type='C',
            molecule_type=type,
            min_pair_cov=min_pair_cov,
            min_base_cov=min_base_cov,
            min_base_rate=min_base_rate,
            fdr_threshold=fdr_threshold,
            logger=logger
        )
    
    logger.info("=" * 60)
    logger.info("Pipeline completed successfully!")
    logger.info("=" * 60)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Double-ended sequencing data processing pipeline"
    )
    parser.add_argument("--r1", required=True, type=Path, help="R1 FASTQ file path")
    parser.add_argument("--r2", required=True, type=Path, help="R2 FASTQ file path")
    parser.add_argument("--genome-fa", required=True, type=Path, help="Genome FASTA file path")
    parser.add_argument("--trans-fa", required=True, type=Path, help="Transcriptome FASTA file path")
    parser.add_argument("--gtf", required=True, type=Path, help="GTF annotation file path")
    parser.add_argument("--tooldir", required=True, type=Path, help="Tool directory (project root)")
    parser.add_argument("--outputdir", required=True, type=Path, help="Output directory")
    parser.add_argument("--prefix", required=True, type=str, help="Output file prefix")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads (default: 8)")
    parser.add_argument("--type", required=True, type=str, choices=['DNA', 'RNA'], 
                       help="Molecule type: DNA or RNA")
    parser.add_argument("--dryrun", action="store_true",
                       help="Dry run mode: print what would be done without executing")
    
    # Filter参数
    parser.add_argument("--filter-target-bases", nargs='+', default=None,
                       help="Target bases for filtering (default: T G for DNA, A C for RNA)")
    parser.add_argument("--filter-cutoffs", nargs='+', type=int, default=None,
                       help="Cutoffs for each target base (default: 3 for all)")
    
    # Call site参数
    parser.add_argument("--min-pair-cov", type=int, default=10,
                       help="Minimum pair coverage for site calling (default: 10)")
    parser.add_argument("--min-base-cov", type=int, default=2,
                       help="Minimum base coverage for site calling (default: 2)")
    parser.add_argument("--min-base-rate", type=float, default=0.1,
                       help="Minimum base rate for site calling (default: 0.1)")
    parser.add_argument("--fdr-threshold", type=float, default=0.05,
                       help="FDR threshold for site calling (default: 0.05)")
    
    args = parser.parse_args()
    
    # 设置日志
    log_dir = args.tooldir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_project_logger(log_dir, f"{args.prefix}.log")
    
    # 运行主流程
    main_pipeline(
        r1_fastq=args.r1,
        r2_fastq=args.r2,
        genome_fa=args.genome_fa,
        trans_fa=args.trans_fa,
        gtf_file=args.gtf,
        tooldir=args.tooldir,
        outputdir=args.outputdir,
        prefix=args.prefix,
        threads=args.threads,
        type=args.type,
        dryrun=args.dryrun,
        logger=logger,
        filter_target_bases=args.filter_target_bases,
        filter_cutoffs=args.filter_cutoffs,
        min_pair_cov=args.min_pair_cov,
        min_base_cov=args.min_base_cov,
        min_base_rate=args.min_base_rate,
        fdr_threshold=args.fdr_threshold
    )

