#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mapping.py

双端测序数据的mapping函数，支持STAR和bowtie2两种工具
"""

import logging
from pathlib import Path
from typing import Optional
from .run_cmd_and_logging import run_cmd, get_logger


def mapping_star(
    r1_fastq: Path,
    r2_fastq: Path,
    index_dir: Path,
    output_prefix: Path,
    type: str = "RNA",
    threads: int = 24,
    mismatch_ratio: float = 0.03,
    logger: Optional[logging.Logger] = None
) -> Path:
    """
    使用STAR进行双端测序数据的mapping
    
    参数:
        r1_fastq: R1 FASTQ文件路径
        r2_fastq: R2 FASTQ文件路径
        index_dir: STAR索引目录路径
        output_prefix: 输出文件前缀（STAR会自动添加后缀，如.bam）
        type: 分子类型，'RNA' 或 'DNA'（默认：'RNA'）
            - RNA: 使用splicing-aware参数
            - DNA: 使用no-splice参数（禁用splicing）
        threads: 线程数（默认：24）
        mismatch_ratio: 错配比例阈值（默认：0.03）
        logger: logger对象（可选）
    
    返回:
        输出的BAM文件路径（Path对象）
    """
    if logger is None:
        logger = get_logger()
    
    # 验证type参数
    type = type.upper()
    if type not in ['RNA', 'DNA']:
        raise ValueError(f"type must be 'RNA' or 'DNA', got '{type}'")
    
    # 确保输出目录存在
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    
    # 根据type选择不同的参数
    if type == 'RNA':
        # RNA模式：使用splicing-aware参数
        cmd_parts = [
            "STAR",
            f"--genomeDir {index_dir}",
            f"--readFilesIn {r1_fastq} {r2_fastq}",
            f"--runThreadN {threads}",
            "--genomeLoad NoSharedMemory",
            f"--outFileNamePrefix {output_prefix}",
            "--outFilterMultimapNmax 20",
            "--alignSJoverhangMin 8",
            "--alignSJDBoverhangMin 1",
            "--outFilterMismatchNmax 2",
            f"--outFilterMismatchNoverReadLmax {mismatch_ratio}",
            "--alignIntronMin 20",
            "--alignIntronMax 1000000",
            "--alignMatesGapMax 1000000",
            "--outSAMunmapped Within",
            "--outFilterType BySJout",
            "--outSAMattributes NH HI AS NM MD",
            "--outSAMtype BAM SortedByCoordinate",
            "--sjdbScore 1",
            "--limitBAMsortRAM 10000000000"
        ]
    else:  # type == 'DNA'
        # DNA模式：使用no-splice参数（禁用splicing）
        cmd_parts = [
            "STAR",
            f"--genomeDir {index_dir}",
            f"--readFilesIn {r1_fastq} {r2_fastq}",
            f"--runThreadN {threads}",
            "--genomeLoad NoSharedMemory",
            f"--outFileNamePrefix {output_prefix}",
            # multimapping
            "--outFilterMultimapNmax 20",
            # ===== NO-SPLICE 核心 =====
            "--alignIntronMax 0",
            "--alignSJoverhangMin 1000000",
            "--alignSJDBoverhangMin 1000000",
            "--outFilterType Normal",
            # mismatch
            "--outFilterMismatchNmax 2",
            f"--outFilterMismatchNoverReadLmax {mismatch_ratio}",
            # PE gap（保留，和 splice 无关）
            "--alignMatesGapMax 1000000",
            # output
            "--outSAMunmapped Within",
            "--outSAMattributes NH HI AS NM MD",
            "--outSAMtype BAM SortedByCoordinate",
            "--limitBAMsortRAM 10000000000"
        ]
    
    # 构建单行命令（shell会自动处理空格）
    cmd = " ".join(cmd_parts)
    
    # 执行STAR命令
    logger.info(f"[STAR] Starting alignment for {r1_fastq.name} and {r2_fastq.name} (type: {type})")
    run_cmd(cmd, logger=logger)
    
    # STAR输出BAM文件路径（SortedByCoordinate会自动生成.bam后缀）
    output_bam = output_prefix.parent / f"{output_prefix.name}Aligned.sortedByCoord.out.bam"
    
    if not output_bam.exists():
        raise FileNotFoundError(f"STAR output BAM file not found: {output_bam}")
    
    logger.info(f"[STAR] Alignment completed. Output BAM: {output_bam}")
    return output_bam


def mapping_bowtie2(
    r1_fastq: Path,
    r2_fastq: Path,
    index_prefix: Path,
    output_sam: Path,
    type: str = "RNA",
    threads: int = 8,
    logger: Optional[logging.Logger] = None
) -> Path:
    """
    使用bowtie2进行双端测序数据的mapping，并自动转换为BAM格式
    
    参数:
        r1_fastq: R1 FASTQ文件路径
        r2_fastq: R2 FASTQ文件路径
        index_prefix: bowtie2索引前缀路径
        output_sam: 输出SAM文件路径（转换后会删除，生成对应的BAM文件）
        type: 分子类型，'RNA' 或 'DNA'（默认：'RNA'）
            - RNA: 正常使用
            - DNA: 会发出警告但仍继续mapping
        threads: 线程数（默认：8）
        logger: logger对象（可选）
    
    返回:
        输出的BAM文件路径（Path对象）
    """
    if logger is None:
        logger = get_logger()
    
    # 验证type参数
    type = type.upper()
    if type not in ['RNA', 'DNA']:
        raise ValueError(f"type must be 'RNA' or 'DNA', got '{type}'")
    
    # 如果是DNA，发出警告但仍继续
    if type == 'DNA':
        logger.warning(
            "[Bowtie2] Warning: Bowtie2 is designed for RNA mapping. "
            "Using DNA type may not be optimal. Proceeding with mapping anyway..."
        )
    
    # 确保输出目录存在
    output_sam.parent.mkdir(parents=True, exist_ok=True)
    
    # 构建bowtie2命令
    cmd_parts = [
        "bowtie2",
        f"-x {index_prefix}",
        f"-1 {r1_fastq}",
        f"-2 {r2_fastq}",
        "--end-to-end",
        "--sensitive",
        "--no-mixed",
        "--no-discordant",
        "-I", "0",
        "-X", "1000",
        f"-p {threads}",
        f"-S {output_sam}"
    ]
    
    # 构建单行命令
    cmd = " ".join(cmd_parts)
    
    # 执行bowtie2命令
    logger.info(f"[Bowtie2] Starting alignment for {r1_fastq.name} and {r2_fastq.name} (type: {type})")
    run_cmd(cmd, logger=logger)
    
    if not output_sam.exists():
        raise FileNotFoundError(f"Bowtie2 output SAM file not found: {output_sam}")
    
    # 将SAM转换为BAM
    output_bam = output_sam.with_suffix('.bam')
    logger.info(f"[Bowtie2] Converting SAM to BAM: {output_sam} -> {output_bam}")
    
    samtools_cmd = f"samtools view -bS {output_sam} -o {output_bam}"
    run_cmd(samtools_cmd, logger=logger)
    
    if not output_bam.exists():
        raise FileNotFoundError(f"BAM file not created: {output_bam}")
    
    # 删除SAM文件
    logger.info(f"[Bowtie2] Removing temporary SAM file: {output_sam}")
    output_sam.unlink()
    
    logger.info(f"[Bowtie2] Alignment completed. Output BAM: {output_bam}")
    return output_bam


def mapping(
    tool: str,
    r1_fastq: Path,
    r2_fastq: Path,
    index_path: Path,
    output_path: Path,
    type: str = "RNA",
    threads: int = 8,
    output_prefix: Optional[Path] = None,
    mismatch_ratio: float = 0.03,
    logger: Optional[logging.Logger] = None
) -> Path:
    """
    双端测序数据的mapping函数，支持STAR和bowtie2
    
    参数:
        tool: 工具名称，'STAR' 或 'bowtie2'
        r1_fastq: R1 FASTQ文件路径
        r2_fastq: R2 FASTQ文件路径
        index_path: 索引路径
            - 对于STAR：索引目录路径
            - 对于bowtie2：索引前缀路径
        output_path: 输出路径
            - 对于STAR：输出文件前缀（会生成.bam文件）
            - 对于bowtie2：输出SAM文件路径（会自动转换为BAM并删除SAM）
        type: 分子类型，'RNA' 或 'DNA'（默认：'RNA'）
            - RNA: 正常使用
            - DNA: 
              * 对于STAR：使用no-splice参数（禁用splicing）
              * 对于bowtie2：会发出警告但仍继续mapping
        threads: 线程数（默认：8）
        output_prefix: 仅用于STAR，输出文件前缀（如果为None，使用output_path）
        mismatch_ratio: 仅用于STAR，错配比例阈值（默认：0.03）
        logger: logger对象（可选）
    
    返回:
        输出的BAM文件路径（Path对象）
    
    示例:
        # STAR mapping (RNA)
        bam = mapping(
            tool='STAR',
            r1_fastq=Path('sample_R1.fq'),
            r2_fastq=Path('sample_R2.fq'),
            index_path=Path('index/STAR_index/'),
            output_path=Path('alignment/sample'),
            type='RNA',
            threads=24,
            mismatch_ratio=0.03
        )
        
        # STAR mapping (DNA)
        bam = mapping(
            tool='STAR',
            r1_fastq=Path('sample_R1.fq'),
            r2_fastq=Path('sample_R2.fq'),
            index_path=Path('index/STAR_index/'),
            output_path=Path('alignment/sample'),
            type='DNA',
            threads=24,
            mismatch_ratio=0.03
        )
        
        # Bowtie2 mapping (RNA)
        bam = mapping(
            tool='bowtie2',
            r1_fastq=Path('sample_R1.fq'),
            r2_fastq=Path('sample_R2.fq'),
            index_path=Path('index/bowtie2_index'),
            output_path=Path('alignment/sample.sam'),
            type='RNA',
            threads=8
        )
    """
    if logger is None:
        logger = get_logger()
    
    tool = tool.upper()
    type = type.upper()
    
    if type not in ['RNA', 'DNA']:
        raise ValueError(f"type must be 'RNA' or 'DNA', got '{type}'")
    
    if tool == 'STAR':
        # 使用output_prefix如果提供，否则使用output_path
        prefix = output_prefix if output_prefix is not None else output_path
        return mapping_star(
            r1_fastq=r1_fastq,
            r2_fastq=r2_fastq,
            index_dir=index_path,
            output_prefix=prefix,
            type=type,
            threads=threads,
            mismatch_ratio=mismatch_ratio,
            logger=logger
        )
    elif tool == 'BOWTIE2':
        return mapping_bowtie2(
            r1_fastq=r1_fastq,
            r2_fastq=r2_fastq,
            index_prefix=index_path,
            output_sam=output_path,
            type=type,
            threads=threads,
            logger=logger
        )
    else:
        raise ValueError(f"Unsupported tool: {tool}. Supported tools: 'STAR', 'bowtie2'")
