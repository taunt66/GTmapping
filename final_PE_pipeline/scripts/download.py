#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
download.py

下载参考基因组、转录组和注释文件，并初始化项目目录结构
"""

import argparse
import gzip
import logging
import shutil
from pathlib import Path
from typing import Optional
import urllib.request
import urllib.error

# 导入日志和命令运行函数
try:
    from .run_cmd_and_logging import run_cmd, get_logger, setup_project_logger
except ImportError:
    from run_cmd_and_logging import run_cmd, get_logger, setup_project_logger


def download_file(url: str, output_path: Path, logger: Optional[logging.Logger] = None) -> Path:
    """
    下载文件（支持断点续传）
    
    参数:
        url: 下载URL
        output_path: 输出文件路径
        logger: logger对象（可选）
    
    返回:
        输出文件路径（Path对象）
    """
    if logger is None:
        logger = get_logger()
    
    # 确保输出目录存在
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # 检查文件是否已存在
    if output_path.exists():
        logger.info(f"File already exists: {output_path}, skipping download")
        return output_path
    
    logger.info(f"Downloading {url} to {output_path}")
    
    try:
        # 使用urllib下载
        def reporthook(blocknum, blocksize, totalsize):
            if totalsize > 0:
                percent = min(100, (blocknum * blocksize * 100) // totalsize)
                if blocknum % 100 == 0:  # 每100个块打印一次
                    logger.info(f"Download progress: {percent}%")
        
        urllib.request.urlretrieve(url, output_path, reporthook=reporthook)
        logger.info(f"Download completed: {output_path}")
        
    except urllib.error.URLError as e:
        logger.error(f"Failed to download {url}: {e}")
        raise
    
    return output_path


def decompress_file(compressed_path: Path, output_path: Optional[Path] = None,
                   logger: Optional[logging.Logger] = None) -> Path:
    """
    解压.gz文件
    
    参数:
        compressed_path: 压缩文件路径
        output_path: 输出文件路径（如果为None，则自动生成）
        logger: logger对象（可选）
    
    返回:
        解压后的文件路径（Path对象）
    """
    if logger is None:
        logger = get_logger()
    
    if output_path is None:
        # 自动生成输出文件名（去掉.gz后缀）
        output_path = compressed_path.with_suffix('')
    
    # 检查输出文件是否已存在
    if output_path.exists():
        logger.info(f"Decompressed file already exists: {output_path}, skipping decompression")
        return output_path
    
    logger.info(f"Decompressing {compressed_path} to {output_path}")
    
    try:
        with gzip.open(compressed_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        logger.info(f"Decompression completed: {output_path}")
    except Exception as e:
        logger.error(f"Failed to decompress {compressed_path}: {e}")
        raise
    
    return output_path


def initialize_project_directories(tooldir: Path, logger: Optional[logging.Logger] = None):
    """
    初始化项目目录结构
    
    参数:
        tooldir: 项目根目录
        logger: logger对象（可选）
    """
    if logger is None:
        logger = get_logger()
    
    logger.info(f"Initializing project directories in {tooldir}")
    
    # 定义目录结构
    directories = {
        'reference': '参考序列目录（基因组、转录组FASTA）',
        'annotations': '注释文件目录（GTF等）',
        'index': '索引文件目录（STAR、bowtie2索引）',
        'logs': '日志文件目录',
        'fastq_data': 'FASTQ数据目录（可选）',
        'alignment': '比对结果目录（可选）'
    }
    
    for dir_name, description in directories.items():
        dir_path = tooldir / dir_name
        dir_path.mkdir(parents=True, exist_ok=True)
        logger.info(f"Created directory: {dir_path} ({description})")
    
    logger.info("Project directories initialized successfully")


def download_ecoli_references(tooldir: Path, logger: Optional[logging.Logger] = None):
    """
    下载E.coli K12 MG1655参考文件
    
    参数:
        tooldir: 项目根目录
        logger: logger对象（可选）
    
    返回:
        (genome_fa, gtf_file) 元组
    """
    if logger is None:
        logger = get_logger()
    
    logger.info("=" * 60)
    logger.info("Downloading E.coli K12 MG1655 references")
    logger.info("=" * 60)
    
    ref_dir = tooldir / "reference"
    anno_dir = tooldir / "annotations"
    
    # 下载URL
    genome_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
    gtf_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gtf.gz"
    
    # 下载文件
    genome_gz = ref_dir / "GCF_000005845.2_ASM584v2_genomic.fna.gz"
    gtf_gz = anno_dir / "GCF_000005845.2_ASM584v2_genomic.gtf.gz"
    
    download_file(genome_url, genome_gz, logger)
    download_file(gtf_url, gtf_gz, logger)
    
    # 解压文件
    genome_fa = decompress_file(genome_gz, logger=logger)
    gtf_file = decompress_file(gtf_gz, logger=logger)
    
    logger.info(f"E.coli genome: {genome_fa}")
    logger.info(f"E.coli GTF: {gtf_file}")
    
    return genome_fa, gtf_file


def download_human_references(tooldir: Path, logger: Optional[logging.Logger] = None):
    """
    下载人类参考文件（基因组、转录组、注释）
    
    参数:
        tooldir: 项目根目录
        logger: logger对象（可选）
    
    返回:
        (genome_fa, trans_fa, gtf_file) 元组
    """
    if logger is None:
        logger = get_logger()
    
    logger.info("=" * 60)
    logger.info("Downloading Human (GRCh38) references")
    logger.info("=" * 60)
    
    ref_dir = tooldir / "reference"
    anno_dir = tooldir / "annotations"
    
    # 下载URL
    genome_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.p14.genome.fa.gz"
    trans_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz"
    gtf_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz"
    
    # 下载文件
    genome_gz = ref_dir / "GRCh38.p14.genome.fa.gz"
    trans_gz = ref_dir / "gencode.v49.transcripts.fa.gz"
    gtf_gz = anno_dir / "gencode.v49.basic.annotation.gtf.gz"
    
    download_file(genome_url, genome_gz, logger)
    download_file(trans_url, trans_gz, logger)
    download_file(gtf_url, gtf_gz, logger)
    
    # 解压文件
    genome_fa = decompress_file(genome_gz, logger=logger)
    trans_fa = decompress_file(trans_gz, logger=logger)
    gtf_file = decompress_file(gtf_gz, logger=logger)
    
    logger.info(f"Human genome: {genome_fa}")
    logger.info(f"Human transcriptome: {trans_fa}")
    logger.info(f"Human GTF: {gtf_file}")
    
    return genome_fa, trans_fa, gtf_file


def main():
    """
    主函数：下载参考文件并初始化项目目录
    """
    parser = argparse.ArgumentParser(
        description="Download reference files and initialize project directories"
    )
    parser.add_argument("--tooldir", required=True, type=Path,
                       help="Project root directory (tool directory)")
    parser.add_argument("--species", required=True, choices=['ecoli', 'human', 'both'],
                       help="Species to download: ecoli, human, or both")
    parser.add_argument("--skip-init", action="store_true",
                       help="Skip directory initialization (if directories already exist)")
    
    args = parser.parse_args()
    
    # 设置日志
    log_dir = args.tooldir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_project_logger(log_dir, "download.log")
    
    logger.info("=" * 60)
    logger.info("Starting download and initialization")
    logger.info("=" * 60)
    logger.info(f"Project directory: {args.tooldir}")
    logger.info(f"Species: {args.species}")
    
    # 初始化项目目录
    if not args.skip_init:
        initialize_project_directories(args.tooldir, logger)
    else:
        logger.info("Skipping directory initialization (--skip-init)")
    
    # 下载参考文件
    if args.species in ['ecoli', 'both']:
        download_ecoli_references(args.tooldir, logger)
    
    if args.species in ['human', 'both']:
        download_human_references(args.tooldir, logger)
    
    logger.info("=" * 60)
    logger.info("Download and initialization completed successfully!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()

