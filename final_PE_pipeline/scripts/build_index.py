#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_index.py

使用不同工具（STAR、bowtie2）构建基因组/转录组索引
"""

from pathlib import Path
from typing import Optional
import math
import logging
from Bio import SeqIO

# 导入日志和命令运行函数
try:
    from .run_cmd_and_logging import run_cmd, get_logger
except ImportError:
    from run_cmd_and_logging import run_cmd, get_logger


def calculate_genome_length(fasta_path: Path) -> int:
    """
    计算FASTA文件中的总序列长度
    
    参数:
        fasta_path: FASTA文件路径
    
    返回:
        总序列长度（bp）
    """
    total_length = 0
    with open(fasta_path, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            total_length += len(record.seq)
    return total_length


def build_star_index(genome_fasta: Path, 
                     output_dir: Path, threads: int = 8,
                     gtf_file: Optional[Path] = None,
                     sjdb_overhang: int = 100,
                     logger: Optional[logging.Logger] = None) -> Path:
    """
    使用STAR构建基因组索引
    
    参数:
        genome_fasta: 基因组FASTA文件路径
        output_dir: 输出索引目录
        threads: 线程数（默认：8）
        gtf_file: GTF注释文件路径（可选，如果提供则用于splice junction数据库）
        sjdb_overhang: splice junction数据库overhang长度（默认：100，仅在提供GTF时使用）
        logger: logger对象（可选）
    
    返回:
        输出索引目录路径（Path对象）
    """
    if logger is None:
        logger = get_logger()
    
    # 确保输出目录存在
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 计算基因组长度
    logger.info(f"Calculating genome length from {genome_fasta}")
    genome_length = calculate_genome_length(genome_fasta)
    logger.info(f"Genome length: {genome_length:,} bp")
    
    # 计算genomeSAindexNbases
    # min(14, log2(GenomeLength)/2 - 1)
    log2_genome = math.log2(genome_length)
    sa_index_nbases = min(14, int(log2_genome / 2 - 1))
    logger.info(f"genomeSAindexNbases: {sa_index_nbases}")
    
    # 构建STAR命令
    cmd_parts = [
        f"STAR",
        f"--runThreadN {threads}",
        f"--runMode genomeGenerate",
        f"--genomeDir {output_dir}",
        f"--genomeFastaFiles {genome_fasta}",
        f"--genomeSAindexNbases {sa_index_nbases}"
    ]
    
    # 如果提供了GTF文件，添加splice junction相关参数
    if gtf_file is not None:
        cmd_parts.append(f"--sjdbGTFfile {gtf_file}")
        cmd_parts.append(f"--sjdbOverhang {sjdb_overhang}")
        logger.info(f"Using GTF file for splice junction database: {gtf_file}")
    else:
        logger.info("Building STAR index without GTF file (no splice junction database)")
    
    cmd = " ".join(cmd_parts)
    
    # 执行命令
    logger.info(f"Building STAR index in {output_dir}")
    run_cmd(cmd, logger=logger, check=True)
    
    logger.info(f"STAR index built successfully in {output_dir}")
    return output_dir


def build_bowtie2_index(genome_fasta: Path, output_prefix: Path,
                        threads: int = 1,
                        logger: Optional[logging.Logger] = None) -> Path:
    """
    使用bowtie2构建基因组索引
    
    参数:
        genome_fasta: 基因组FASTA文件路径
        output_prefix: 输出索引文件前缀（不含扩展名）
        threads: 线程数（默认：1）
        logger: logger对象（可选）
    
    返回:
        输出索引文件前缀路径（Path对象）
    
    注意:
        bowtie2-build会生成多个文件（.1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, .rev.2.bt2）
        返回的是前缀路径，实际文件会在此基础上添加扩展名
    """
    if logger is None:
        logger = get_logger()
    
    # 确保输出目录存在
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    
    # 构建bowtie2-build命令
    cmd = (
        f"bowtie2-build "
        f"--threads {threads} "
        f"{genome_fasta} "
        f"{output_prefix}"
    )
    
    # 执行命令
    logger.info(f"Building bowtie2 index with prefix {output_prefix}")
    run_cmd(cmd, logger=logger, check=True)
    
    logger.info(f"bowtie2 index built successfully with prefix {output_prefix}")
    return output_prefix


def build_index(tool: str, genome_fasta: Path, output_path: Path,
                threads: int = 8, gtf_file: Optional[Path] = None,
                sjdb_overhang: int = 100,
                logger: Optional[logging.Logger] = None) -> Path:
    """
    使用指定工具构建索引（统一接口）
    
    参数:
        tool: 工具名称，'STAR' 或 'bowtie2'
        genome_fasta: 基因组FASTA文件路径
        output_path: 输出路径
            - STAR: 索引目录路径
            - bowtie2: 索引文件前缀路径
        threads: 线程数（默认：8）
        gtf_file: GTF注释文件路径（可选，如果提供则用于STAR的splice junction数据库）
        sjdb_overhang: STAR的splice junction overhang长度（默认：100，仅在提供GTF时使用）
        logger: logger对象（可选）
    
    返回:
        输出路径（Path对象）
    
    示例:
        # STAR索引（带GTF）
        build_index('STAR', Path('genome.fa'), Path('star_index/'), 
                   threads=16, gtf_file=Path('annotation.gtf'))
        
        # STAR索引（不带GTF）
        build_index('STAR', Path('genome.fa'), Path('star_index/'), 
                   threads=16)
        
        # bowtie2索引
        build_index('bowtie2', Path('genome.fa'), Path('bowtie2_index/genome'),
                   threads=30)
    """
    if logger is None:
        logger = get_logger()
    
    tool = tool.upper()
    
    if tool == 'STAR':
        return build_star_index(
            genome_fasta=genome_fasta,
            output_dir=output_path,
            threads=threads,
            gtf_file=gtf_file,
            sjdb_overhang=sjdb_overhang,
            logger=logger
        )
    elif tool == 'BOWTIE2':
        return build_bowtie2_index(
            genome_fasta=genome_fasta,
            output_prefix=output_path,
            threads=threads,
            logger=logger
        )
    else:
        raise ValueError(f"Unsupported tool: {tool}. Supported tools: 'STAR', 'bowtie2'")


if __name__ == "__main__":
    import sys
    from pathlib import Path
    
    # 示例用法
    if len(sys.argv) < 4:
        print("Usage: python build_index.py <tool> <genome_fasta> <output_path> [threads] [gtf_file]")
        print("  tool: STAR or bowtie2")
        print("  genome_fasta: Path to genome FASTA file")
        print("  output_path: Output directory (STAR) or prefix (bowtie2)")
        print("  threads: Number of threads (default: 8)")
        print("  gtf_file: GTF file (required for STAR)")
        sys.exit(1)
    
    tool = sys.argv[1]
    genome_fasta = Path(sys.argv[2])
    output_path = Path(sys.argv[3])
    threads = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    gtf_file = Path(sys.argv[5]) if len(sys.argv) > 5 else None
    
    # 设置logger（如果需要）
    try:
        from run_cmd_and_logging import setup_project_logger
        log_dir = Path("./logs")
        logger = setup_project_logger(log_dir, "build_index.log")
    except ImportError:
        logger = None
    
    build_index(tool, genome_fasta, output_path, threads, gtf_file, logger=logger)

