#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spike_in_analysis.py

分析spike-in序列的甲基化水平

流程：
1. 读取spike-in FASTA文件和index文件
2. 检查并构建bowtie2索引（如果不存在）
3. 根据type转换FASTQ文件（与run.py一致）
4. 使用bowtie2进行mapping
5. 进行pileup分析
6. 统计目标位点的甲基化水平
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import defaultdict

# 处理导入路径
_script_dir = Path(__file__).parent
_scripts_dir = _script_dir / "scripts"

if str(_scripts_dir) not in sys.path:
    sys.path.insert(0, str(_scripts_dir))
if str(_script_dir.parent) not in sys.path:
    sys.path.insert(0, str(_script_dir.parent))

# 尝试导入
try:
    from final_PE_pipeline.scripts.build_index import build_index
    from final_PE_pipeline.scripts.convert_fastq_fasta import change_fastq, change_fasta
    from final_PE_pipeline.scripts.mapping import mapping
    from final_PE_pipeline.scripts.pileup_PE import get_genome_wide_pileup_paired
    from final_PE_pipeline.scripts.run_cmd_and_logging import setup_project_logger, get_logger, run_cmd
except ImportError:
    from scripts.build_index import build_index
    from scripts.convert_fastq_fasta import change_fastq, change_fasta
    from scripts.mapping import mapping
    from scripts.pileup_PE import get_genome_wide_pileup_paired
    from scripts.run_cmd_and_logging import setup_project_logger, get_logger, run_cmd


def parse_index_file(index_file: Path) -> Dict[str, List[int]]:
    """
    解析index文件，返回序列ID到目标碱基位置的映射
    
    参数:
        index_file: index文件路径
        
    返回:
        Dict[序列ID, 目标碱基位置列表（1-based）]
    """
    index_dict = {}
    current_id = None
    
    with open(index_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                # 序列ID行
                current_id = line[1:]  # 去掉>符号
            else:
                # 位置行
                if current_id is None:
                    raise ValueError(f"Index file format error: position line without sequence ID")
                
                # 解析位置（可能用_分隔多个位置）
                positions = []
                for pos_str in line.split('_'):
                    try:
                        pos = int(pos_str)
                        positions.append(pos)
                    except ValueError:
                        raise ValueError(f"Invalid position in index file: {pos_str}")
                
                index_dict[current_id] = positions
                current_id = None  # 重置，准备下一个序列
    
    return index_dict


def analyze_spikein_methylation(
    spikein_fa: Path,
    index_file: Path,
    r1_fastq: Path,
    r2_fastq: Path,
    tooldir: Path,
    outputdir: Path,
    prefix: str,
    type: str,
    threads: int = 8,
    logger: Optional[logging.Logger] = None
) -> Path:
    """
    分析spike-in序列的甲基化水平
    
    参数:
        spikein_fa: spike-in FASTA文件路径
        index_file: index文件路径（包含目标碱基位置）
        r1_fastq: R1 FASTQ文件路径
        r2_fastq: R2 FASTQ文件路径
        tooldir: 工具目录（项目根目录）
        outputdir: 输出目录
        prefix: 输出文件前缀
        type: 分子类型，'DNA' 或 'RNA'
        threads: 线程数（默认：8）
        logger: logger对象（可选）
    
    返回:
        甲基化水平统计结果文件路径
    """
    if logger is None:
        logger = get_logger()
    
    type = type.upper()
    if type not in ['DNA', 'RNA']:
        raise ValueError(f"type must be 'DNA' or 'RNA', got '{type}'")
    
    logger.info("=" * 60)
    logger.info("Spike-in Methylation Analysis")
    logger.info("=" * 60)
    logger.info(f"Spike-in FASTA: {spikein_fa}")
    logger.info(f"Index file: {index_file}")
    logger.info(f"R1 FASTQ: {r1_fastq}")
    logger.info(f"R2 FASTQ: {r2_fastq}")
    logger.info(f"Type: {type}")
    logger.info(f"Threads: {threads}")
    logger.info("=" * 60)
    
    # 创建输出目录
    outputdir.mkdir(parents=True, exist_ok=True)
    index_dir = tooldir / "index"
    index_dir.mkdir(parents=True, exist_ok=True)
    fastq_dir = outputdir / "fastq"
    fastq_dir.mkdir(parents=True, exist_ok=True)
    alignment_dir = outputdir / "alignment"
    alignment_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: 解析index文件
    logger.info("=" * 60)
    logger.info("Step 1: Parsing index file")
    logger.info("=" * 60)
    index_dict = parse_index_file(index_file)
    logger.info(f"Found {len(index_dict)} spike-in sequences with target positions")
    for seq_id, positions in index_dict.items():
        logger.info(f"  {seq_id}: positions {positions}")
    
    # Step 2: 转换spike-in FASTA文件（根据type）
    logger.info("=" * 60)
    logger.info("Step 2: Converting spike-in FASTA file")
    logger.info("=" * 60)
    
    # 根据type确定转换规则（与FASTQ转换规则一致）
    # DNA: R1转换是T->C, G->A，所以FASTA也应该做T->C, G->A转换
    # RNA: R1转换是A->G, C->T，所以FASTA也应该做A->G, C->T转换
    if type == 'DNA':
        # DNA: FASTA进行T->C, G->A转换（与R1转换一致）
        spikein_convert_list = ['T:C', 'G:A']
        spikein_convert_label = "TC_GA"
    else:  # RNA
        # RNA: FASTA进行A->G, C->T转换（与R1转换一致）
        spikein_convert_list = ['A:G', 'C:T']
        spikein_convert_label = "AG_CT"
    
    # 定义转换后的FASTA文件路径
    # 注意：spike-in FASTA可以在任意位置，转换后的文件会保存在tooldir/reference/下
    ref_dir = tooldir / "reference"
    ref_dir.mkdir(parents=True, exist_ok=True)
    
    # 使用输入文件的stem（不包含路径）来构建输出文件名
    spikein_converted_fa = ref_dir / f"{spikein_fa.stem}_{spikein_convert_label}convert.fa"
    
    logger.info(f"Input spike-in FASTA: {spikein_fa}")
    logger.info(f"Output converted FASTA: {spikein_converted_fa}")
    logger.info(f"Conversion rules: {', '.join(spikein_convert_list)}")
    
    if not spikein_converted_fa.exists():
        logger.info(f"Converting spike-in FASTA ({type}): {', '.join(spikein_convert_list)}...")
        change_fasta(
            input_fasta=spikein_fa,
            output_fasta=spikein_converted_fa,
            convert_list=spikein_convert_list
        )
        logger.info(f"Converted spike-in FASTA: {spikein_converted_fa}")
    else:
        logger.info(f"Converted spike-in FASTA already exists: {spikein_converted_fa}")
    
    # Step 3: 检查并构建bowtie2索引（使用转换后的FASTA）
    logger.info("=" * 60)
    logger.info("Step 3: Building bowtie2 index")
    logger.info("=" * 60)
    spikein_index_prefix = index_dir / f"{spikein_fa.stem}_{spikein_convert_label}_bowtie2_index"
    
    logger.info(f"Index prefix: {spikein_index_prefix}")
    logger.info(f"Using converted FASTA for indexing: {spikein_converted_fa}")
    
    if not spikein_index_prefix.exists() or not (spikein_index_prefix.with_suffix('.1.bt2')).exists():
        logger.info(f"Building bowtie2 index for converted spike-in sequences...")
        build_index(
            tool='bowtie2',
            genome_fasta=spikein_converted_fa,  # 使用转换后的FASTA
            output_path=spikein_index_prefix,
            threads=threads,
            logger=logger
        )
        logger.info(f"Bowtie2 index: {spikein_index_prefix}")
    else:
        logger.info(f"Bowtie2 index already exists: {spikein_index_prefix}")
    
    # Step 4: 转换FASTQ文件（根据type）
    logger.info("=" * 60)
    logger.info("Step 4: Converting FASTQ files")
    logger.info("=" * 60)
    
    # 根据type确定转换规则（与run.py一致）
    if type == 'DNA':
        # DNA: R1进行T->C, G->A转换；R2进行A->G, C->T转换
        r1_convert_list = ['T:C', 'G:A']
        r2_convert_list = ['A:G', 'C:T']
        r1_index_label = "T-C_G-A"
        r2_index_label = "A-G_C-T"
    else:  # RNA
        # RNA: R1进行A->G, C->T转换；R2进行T->C, G->A转换
        r1_convert_list = ['A:G', 'C:T']
        r2_convert_list = ['T:C', 'G:A']
        r1_index_label = "A-G_C-T"
        r2_index_label = "T-C_G-A"
    
    # R1转换
    r1_converted_fq = fastq_dir / f"{prefix}_r1_converted.fq.gz"
    r1_index_file = fastq_dir / f"{prefix}_r1_index_{r1_index_label}.txt"
    
    logger.info(f"Converting R1 FASTQ ({type}): {', '.join(r1_convert_list)}...")
    r1_converted, r1_idx = change_fastq(
        input_fastq=r1_fastq,
        output_dir=fastq_dir,
        convert_list=r1_convert_list,
        output_fastq=r1_converted_fq
    )
    if r1_idx != r1_index_file:
        if r1_idx.exists():
            r1_idx.rename(r1_index_file)
        else:
            logger.warning(f"Expected index file {r1_index_file}, but got {r1_idx}")
            r1_index_file = r1_idx
    logger.info(f"R1 converted FASTQ: {r1_converted_fq}")
    logger.info(f"R1 index file: {r1_index_file}")
    
    # R2转换
    r2_converted_fq = fastq_dir / f"{prefix}_r2_converted.fq.gz"
    r2_index_file = fastq_dir / f"{prefix}_r2_index_{r2_index_label}.txt"
    
    logger.info(f"Converting R2 FASTQ ({type}): {', '.join(r2_convert_list)}...")
    r2_converted, r2_idx = change_fastq(
        input_fastq=r2_fastq,
        output_dir=fastq_dir,
        convert_list=r2_convert_list,
        output_fastq=r2_converted_fq
    )
    if r2_idx != r2_index_file:
        if r2_idx.exists():
            r2_idx.rename(r2_index_file)
        else:
            logger.warning(f"Expected index file {r2_index_file}, but got {r2_idx}")
            r2_index_file = r2_idx
    logger.info(f"R2 converted FASTQ: {r2_converted_fq}")
    logger.info(f"R2 index file: {r2_index_file}")
    
    # Step 5: 使用bowtie2进行mapping
    logger.info("=" * 60)
    logger.info("Step 5: Mapping with bowtie2")
    logger.info("=" * 60)
    
    bam_sam = alignment_dir / f"{prefix}.bowtie2_spikein.sam"
    bam_output = alignment_dir / f"{prefix}.bowtie2_spikein.bam"
    
    logger.info(f"Mapping to spike-in sequences...")
    bam_output = mapping(
        tool='bowtie2',
        r1_fastq=r1_converted_fq,
        r2_fastq=r2_converted_fq,
        index_path=spikein_index_prefix,
        output_path=bam_sam,
        type=type,
        threads=threads,
        logger=logger
    )
    logger.info(f"Mapped BAM (sorted & indexed by mapping function): {bam_output}")
    
    # Step 6: 进行pileup分析
    logger.info("=" * 60)
    logger.info("Step 6: Generating pileup")
    logger.info("=" * 60)
    
    # 根据type确定主要read和允许的ref bases
    if type == 'DNA':
        primary_read = 'R2'  # DNA按R2计数
        allowed_ref_bases = {'A', 'C'}  # DNA统计A和C
    else:  # RNA
        primary_read = 'R1'  # RNA按R1计数
        allowed_ref_bases = {'A', 'C'}  # RNA统计A和C
    
    pileup_file = outputdir / f"{prefix}.pileup"
    
    logger.info(f"Generating pileup (primary_read={primary_read}, allowed_ref_bases={allowed_ref_bases})...")
    get_genome_wide_pileup_paired(
        bam_path=bam_output,
        fasta_path=spikein_converted_fa,  # 使用转换后的spike-in FASTA（与mapping一致）
        outputdir=outputdir,
        thread=threads,
        prefix=prefix,
        tmp_dir=None,
        chunk_size=500000,
        max_depth=100000,
        allowed_ref_bases=allowed_ref_bases,
        primary_read=primary_read,
        logger=logger
    )
    logger.info(f"Pileup file: {pileup_file}")
    
    # Step 7: 统计目标位点的甲基化水平
    logger.info("=" * 60)
    logger.info("Step 7: Analyzing methylation levels at target positions")
    logger.info("=" * 60)
    
    result_file = outputdir / f"{prefix}.spikein_methylation.txt"
    
    # 读取pileup文件并统计
    methylation_stats = analyze_pileup_for_positions(
        pileup_file=pileup_file,
        index_dict=index_dict,
        spikein_fa=spikein_fa,
        result_file=result_file,
        logger=logger
    )
    
    logger.info("=" * 60)
    logger.info("Spike-in methylation analysis completed!")
    logger.info(f"Results saved to: {result_file}")
    logger.info("=" * 60)
    
    return result_file


def analyze_pileup_for_positions(
    pileup_file: Path,
    index_dict: Dict[str, List[int]],
    spikein_fa: Path,
    result_file: Path,
    logger: Optional[logging.Logger] = None
) -> Dict[str, Dict[int, Dict[str, float]]]:
    """
    从pileup文件中统计目标位点的甲基化水平
    
    参数:
        pileup_file: pileup文件路径
        index_dict: 序列ID到目标碱基位置的映射
        spikein_fa: spike-in FASTA文件路径（用于读取参考碱基）
        result_file: 输出结果文件路径
        logger: logger对象（可选）
    
    返回:
        甲基化统计结果字典
    """
    if logger is None:
        logger = get_logger()
    
    # 读取FASTA文件，获取参考序列
    from Bio import SeqIO
    ref_sequences = {}
    with open(spikein_fa, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            ref_sequences[record.id] = str(record.seq).upper()
    
    # 读取pileup文件
    # pileup格式：chrom  pos  ref_base  strand  bases
    # 注意：pos是1-based
    
    # 存储每个位点的统计信息
    # stats[seq_id][pos] = {'A': count, 'G': count, 'C': count, 'T': count, 'total': count}
    stats = defaultdict(lambda: defaultdict(lambda: {'A': 0, 'G': 0, 'C': 0, 'T': 0, 'total': 0}))
    
    logger.info(f"Reading pileup file: {pileup_file}")
    with open(pileup_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            fields = line.split('\t')
            if len(fields) < 5:
                continue
            
            chrom = fields[0]
            pos = int(fields[1])
            ref_base = fields[2]
            strand = fields[3]
            bases = fields[4]
            
            # 检查这个位点是否在index_dict中
            if chrom in index_dict and pos in index_dict[chrom]:
                # 统计碱基计数
                base_counts = {'A': 0, 'G': 0, 'C': 0, 'T': 0}
                total = 0
                
                for base in bases:
                    base_upper = base.upper()
                    if base_upper in base_counts:
                        base_counts[base_upper] += 1
                        total += 1
                
                # 更新统计
                stats[chrom][pos]['A'] += base_counts['A']
                stats[chrom][pos]['G'] += base_counts['G']
                stats[chrom][pos]['C'] += base_counts['C']
                stats[chrom][pos]['T'] += base_counts['T']
                stats[chrom][pos]['total'] += total
    
    # 写入结果文件
    logger.info(f"Writing results to: {result_file}")
    with open(result_file, 'w') as f:
        # 写入表头
        f.write("Sequence_ID\tPosition\tRef_Base\tA_count\tG_count\tC_count\tT_count\tTotal_count\tA_rate\tG_rate\tC_rate\tT_rate\n")
        
        # 写入每个位点的统计
        for seq_id in sorted(index_dict.keys()):
            positions = index_dict[seq_id]
            # 获取参考序列
            ref_seq = ref_sequences.get(seq_id, '')
            
            for pos in sorted(positions):
                if seq_id in stats and pos in stats[seq_id]:
                    stat = stats[seq_id][pos]
                    total = stat['total']
                    
                    # 计算各碱基的比例
                    if total > 0:
                        a_rate = stat['A'] / total
                        g_rate = stat['G'] / total
                        c_rate = stat['C'] / total
                        t_rate = stat['T'] / total
                    else:
                        a_rate = g_rate = c_rate = t_rate = 0.0
                    
                    # 获取参考碱基（从FASTA文件中读取，1-based位置）
                    if ref_seq and 1 <= pos <= len(ref_seq):
                        ref_base = ref_seq[pos - 1]  # 转换为0-based索引
                    else:
                        ref_base = 'N'
                    
                    f.write(f"{seq_id}\t{pos}\t{ref_base}\t"
                           f"{stat['A']}\t{stat['G']}\t{stat['C']}\t{stat['T']}\t{total}\t"
                           f"{a_rate:.4f}\t{g_rate:.4f}\t{c_rate:.4f}\t{t_rate:.4f}\n")
                else:
                    # 如果该位点没有覆盖，写入0
                    # 获取参考碱基
                    if ref_seq and 1 <= pos <= len(ref_seq):
                        ref_base = ref_seq[pos - 1]
                    else:
                        ref_base = 'N'
                    
                    f.write(f"{seq_id}\t{pos}\t{ref_base}\t0\t0\t0\t0\t0\t0.0000\t0.0000\t0.0000\t0.0000\n")
    
    return stats


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analyze spike-in methylation levels"
    )
    parser.add_argument("--spikein-fa", required=True, type=Path,
                       help="Spike-in FASTA file path")
    parser.add_argument("--index", required=True, type=Path,
                       help="Index file path (contains target base positions)")
    parser.add_argument("--r1", required=True, type=Path,
                       help="R1 FASTQ file path")
    parser.add_argument("--r2", required=True, type=Path,
                       help="R2 FASTQ file path")
    parser.add_argument("--tooldir", required=True, type=Path,
                       help="Tool directory (project root)")
    parser.add_argument("--outputdir", required=True, type=Path,
                       help="Output directory")
    parser.add_argument("--prefix", required=True, type=str,
                       help="Output file prefix")
    parser.add_argument("--type", required=True, type=str, choices=['DNA', 'RNA'],
                       help="Molecule type: DNA or RNA")
    parser.add_argument("--threads", type=int, default=8,
                       help="Number of threads (default: 8)")
    
    args = parser.parse_args()
    
    # 设置日志
    log_dir = args.tooldir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_project_logger(log_dir, f"{args.prefix}_spikein.log")
    
    # 运行分析
    analyze_spikein_methylation(
        spikein_fa=args.spikein_fa,
        index_file=args.index,
        r1_fastq=args.r1,
        r2_fastq=args.r2,
        tooldir=args.tooldir,
        outputdir=args.outputdir,
        prefix=args.prefix,
        type=args.type,
        threads=args.threads,
        logger=logger
    )

