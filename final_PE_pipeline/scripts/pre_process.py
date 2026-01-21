#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pre_process.py

双端测序数据预处理脚本

功能：
- 自动识别输入输出是否为压缩格式（根据.gz后缀）
- 根据type参数（DNA或RNA）进行不同处理：
  - RNA: 直接pass（不做处理）
  - DNA: 
    - R1: 如果开头是9个G，则取第10个碱基之后的内容
    - R2: 使用cutadapt -a {C}9 去除3'端的9个C
"""

import argparse
import gzip
import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

# 导入日志和命令运行函数
try:
    from .run_cmd_and_logging import run_cmd, get_logger, setup_project_logger
except ImportError:
    from run_cmd_and_logging import run_cmd, get_logger, setup_project_logger


def check_tool(name: str) -> str:
    """
    检查工具是否在PATH中
    
    参数:
        name: 工具名称
    
    返回:
        工具路径
    
    异常:
        RuntimeError: 如果工具未找到
    """
    tool_path = shutil.which(name)
    if tool_path is None:
        raise RuntimeError(f"Tool {name} not found in PATH. Please install it first.")
    return tool_path


def open_read(path: Path):
    """
    根据文件扩展名自动选择打开方式（支持.gz）
    
    参数:
        path: 文件路径
    
    返回:
        文件句柄
    """
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def open_write(path: Path):
    """
    根据文件扩展名自动选择写入方式（支持.gz）
    
    参数:
        path: 文件路径
    
    返回:
        文件句柄
    """
    # 确保输出目录存在
    path.parent.mkdir(parents=True, exist_ok=True)
    
    if str(path).endswith('.gz'):
        return gzip.open(path, 'wt')
    return open(path, 'w')


def process_r1_dna(r1_in: Path, r1_out: Path, logger: Optional[logging.Logger] = None):
    """
    处理DNA类型的R1文件：如果开头是9个G，则取第10个碱基之后的内容
    
    参数:
        r1_in: 输入R1 FASTQ文件路径
        r1_out: 输出R1 FASTQ文件路径
        logger: logger对象（可选）
    """
    if logger is None:
        logger = get_logger()
    
    logger.info(f"Processing R1 (DNA): checking for 9 G's at the beginning")
    logger.info(f"  Input: {r1_in}")
    logger.info(f"  Output: {r1_out}")
    
    count_total = 0
    count_trimmed = 0
    count_skipped = 0
    
    with open_read(r1_in) as inf, open_write(r1_out) as outf:
        while True:
            # 读取4行（一个FASTQ记录）
            header = inf.readline()
            if not header:
                break
            
            seq = inf.readline().rstrip('\n')
            plus = inf.readline()
            qual = inf.readline().rstrip('\n')
            
            count_total += 1
            
            # 检查开头是否是9个G
            if len(seq) >= 9 and seq[:9] == 'G' * 9:
                # 取第10个碱基之后的内容（索引9开始）
                seq_trimmed = seq[9:]
                qual_trimmed = qual[9:] if len(qual) >= 9 else qual
                
                # 如果trim后序列为空，跳过这条read
                if len(seq_trimmed) == 0:
                    count_skipped += 1
                    continue
                
                # 写入trimmed后的read
                outf.write(header)
                outf.write(seq_trimmed + '\n')
                outf.write(plus)
                outf.write(qual_trimmed + '\n')
                count_trimmed += 1
            else:
                # 开头不是9个G，直接输出原序列
                outf.write(header)
                outf.write(seq + '\n')
                outf.write(plus)
                outf.write(qual + '\n')
    
    logger.info(f"R1 processing completed:")
    logger.info(f"  Total reads: {count_total}")
    logger.info(f"  Trimmed (9 G's found): {count_trimmed}")
    logger.info(f"  Skipped (empty after trim): {count_skipped}")
    logger.info(f"  Unchanged: {count_total - count_trimmed - count_skipped}")


def process_r2_dna(r2_in: Path, r2_out: Path, threads: int = 4,
                  logger: Optional[logging.Logger] = None):
    """
    处理DNA类型的R2文件：使用cutadapt去除3'端的9个C
    
    参数:
        r2_in: 输入R2 FASTQ文件路径
        r2_out: 输出R2 FASTQ文件路径
        threads: 线程数（默认：4）
        logger: logger对象（可选）
    """
    if logger is None:
        logger = get_logger()
    
    logger.info(f"Processing R2 (DNA): removing 9 C's from 3' end using cutadapt")
    logger.info(f"  Input: {r2_in}")
    logger.info(f"  Output: {r2_out}")
    
    cutadapt = check_tool("cutadapt")
    
    # cutadapt命令：-a {C}9 表示去除3'端的9个C
    # -j 指定线程数
    cmd = (
        f"{cutadapt} "
        f"-a C{{9}} "
        f"-j {threads} "
        f"-o {r2_out} "
        f"{r2_in}"
    )
    
    run_cmd(cmd, logger=logger, check=True)
    logger.info(f"R2 processing completed: {r2_out}")


def preprocess_paired_fastq(
    r1_in: Path,
    r2_in: Path,
    r1_out: Path,
    r2_out: Path,
    type: str,
    threads: int = 4,
    logger: Optional[logging.Logger] = None
):
    """
    双端FASTQ预处理主函数
    
    参数:
        r1_in: 输入R1 FASTQ文件路径（自动识别.gz）
        r2_in: 输入R2 FASTQ文件路径（自动识别.gz）
        r1_out: 输出R1 FASTQ文件路径（自动识别.gz）
        r2_out: 输出R2 FASTQ文件路径（自动识别.gz）
        type: 类型，'DNA' 或 'RNA'
        threads: 线程数（默认：4，用于cutadapt）
        logger: logger对象（可选）
    """
    if logger is None:
        logger = get_logger()
    
    type = type.upper()
    if type not in ['DNA', 'RNA']:
        raise ValueError(f"type must be 'DNA' or 'RNA', got '{type}'")
    
    logger.info("=" * 60)
    logger.info("Starting paired-end FASTQ preprocessing")
    logger.info("=" * 60)
    logger.info(f"Type: {type}")
    logger.info(f"Input R1: {r1_in}")
    logger.info(f"Input R2: {r2_in}")
    logger.info(f"Output R1: {r1_out}")
    logger.info(f"Output R2: {r2_out}")
    logger.info(f"Threads: {threads}")
    
    if type == 'RNA':
        logger.info("RNA type: passing through without modification")
        # 直接复制文件
        logger.info("Copying R1...")
        with open_read(r1_in) as inf, open_write(r1_out) as outf:
            shutil.copyfileobj(inf, outf)
        logger.info("Copying R2...")
        with open_read(r2_in) as inf, open_write(r2_out) as outf:
            shutil.copyfileobj(inf, outf)
        logger.info("RNA preprocessing completed (pass-through)")
    
    elif type == 'DNA':
        logger.info("DNA type: processing R1 and R2")
        
        # 处理R1：检查9个G并trim
        process_r1_dna(r1_in, r1_out, logger)
        
        # 处理R2：使用cutadapt去除3'端的9个C
        process_r2_dna(r2_in, r2_out, threads, logger)
        
        logger.info("DNA preprocessing completed")
    
    logger.info("=" * 60)
    logger.info("Preprocessing completed successfully!")
    logger.info("=" * 60)


def main():
    """
    主函数：命令行接口
    """
    parser = argparse.ArgumentParser(
        description="Preprocess paired-end FASTQ files"
    )
    parser.add_argument("-1", "--r1", required=True, type=Path,
                       help="Input R1 FASTQ file (.fq/.fq.gz)")
    parser.add_argument("-2", "--r2", required=True, type=Path,
                       help="Input R2 FASTQ file (.fq/.fq.gz)")
    parser.add_argument("-o1", "--r1-out", required=True, type=Path,
                       help="Output R1 FASTQ file (.fq/.fq.gz)")
    parser.add_argument("-o2", "--r2-out", required=True, type=Path,
                       help="Output R2 FASTQ file (.fq/.fq.gz)")
    parser.add_argument("--type", required=True, choices=['DNA', 'RNA'],
                       help="Type: DNA or RNA")
    parser.add_argument("-t", "--threads", type=int, default=4,
                       help="Number of threads (default: 4)")
    parser.add_argument("--log-dir", type=Path, default=None,
                       help="Log directory (optional)")
    
    args = parser.parse_args()
    
    # 设置日志
    if args.log_dir:
        log_dir = args.log_dir
        log_dir.mkdir(parents=True, exist_ok=True)
        logger = setup_project_logger(log_dir, "preprocess.log")
    else:
        logger = get_logger()
    
    # 运行预处理
    preprocess_paired_fastq(
        r1_in=args.r1,
        r2_in=args.r2,
        r1_out=args.r1_out,
        r2_out=args.r2_out,
        type=args.type,
        threads=args.threads,
        logger=logger
    )


if __name__ == "__main__":
    main()
