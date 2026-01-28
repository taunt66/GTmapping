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
    - Step 1: 使用Python切掉R1 5'端的所有G，根据G数量过滤pair
      * 如果切掉的G数量 >= 7，保留pair并trim掉这些G
      * 如果切掉的G数量 < 7，丢弃pair
    - Step 2: 使用cutadapt去除接头序列
      * R1 3'端: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
      * R2 3'端: CCCCCCC (7个C)
      * 最小长度: 20bp
"""

import argparse
import gzip
import logging
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

# 导入日志和命令运行函数
try:
    from .scripts.run_cmd_and_logging import run_cmd, get_logger, setup_project_logger
except ImportError:
    from scripts.run_cmd_and_logging import run_cmd, get_logger, setup_project_logger


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


def trim_polyg_and_filter_pairs(
    r1_in: Path,
    r2_in: Path,
    r1_out: Path,
    r2_out: Path,
    min_g_count: int = 7,
    logger: Optional[logging.Logger] = None
):
    """
    处理DNA类型的双端FASTQ：切掉R1 5'端的所有G，根据切掉的G数量过滤pair
    
    算法：
    1. 对R1，统计5'端连续G的数量
    2. 如果G数量 >= min_g_count，保留pair，trim掉这些G
    3. 如果G数量 < min_g_count，丢弃pair
    
    参数:
        r1_in: 输入R1 FASTQ文件路径
        r2_in: 输入R2 FASTQ文件路径
        r1_out: 输出R1 FASTQ文件路径（trimmed后）
        r2_out: 输出R2 FASTQ文件路径（与保留的R1配对）
        min_g_count: 最小G数量阈值（默认：7）
        logger: logger对象（可选）
    
    返回:
        (保留的pair数量, 丢弃的pair数量)
    """
    if logger is None:
        logger = get_logger()
    
    logger.info(f"Processing DNA paired-end FASTQ: trimming 5' poly-G from R1")
    logger.info(f"  Input R1: {r1_in}")
    logger.info(f"  Input R2: {r2_in}")
    logger.info(f"  Output R1: {r1_out}")
    logger.info(f"  Output R2: {r2_out}")
    logger.info(f"  Minimum G count threshold: {min_g_count}")
    
    count_total = 0
    count_kept = 0
    count_discarded = 0
    
    with open_read(r1_in) as r1_inf, open_read(r2_in) as r2_inf, \
         open_write(r1_out) as r1_outf, open_write(r2_out) as r2_outf:
        
        while True:
            # 读取R1的4行（一个FASTQ记录）
            r1_header = r1_inf.readline()
            if not r1_header:
                break
            
            r1_seq = r1_inf.readline().rstrip('\n')
            r1_plus = r1_inf.readline()
            r1_qual = r1_inf.readline().rstrip('\n')
            
            # 读取R2的4行（一个FASTQ记录）
            r2_header = r2_inf.readline()
            r2_seq = r2_inf.readline().rstrip('\n')
            r2_plus = r2_inf.readline()
            r2_qual = r2_inf.readline().rstrip('\n')
            
            count_total += 1
            
            # 统计R1 5'端连续G的数量
            g_count = 0
            for base in r1_seq:
                if base == 'G' or base == 'g':
                    g_count += 1
                else:
                    break
            
            # 判断是否保留pair
            if g_count >= min_g_count:
                # 保留pair，trim掉R1的G
                r1_seq_trimmed = r1_seq[g_count:]
                r1_qual_trimmed = r1_qual[g_count:] if len(r1_qual) >= g_count else r1_qual
                
                # 如果trim后序列为空，跳过这条pair
                if len(r1_seq_trimmed) == 0:
                    count_discarded += 1
                    continue
                
                # 写入trimmed后的R1和对应的R2
                r1_outf.write(r1_header)
                r1_outf.write(r1_seq_trimmed + '\n')
                r1_outf.write(r1_plus)
                r1_outf.write(r1_qual_trimmed + '\n')
                
                r2_outf.write(r2_header)
                r2_outf.write(r2_seq + '\n')
                r2_outf.write(r2_plus)
                r2_outf.write(r2_qual + '\n')
                
                count_kept += 1
            else:
                # G数量不足，丢弃pair
                count_discarded += 1
    
    logger.info(f"Poly-G trimming completed:")
    logger.info(f"  Total pairs: {count_total}")
    logger.info(f"  Kept pairs (G count >= {min_g_count}): {count_kept}")
    logger.info(f"  Discarded pairs (G count < {min_g_count}): {count_discarded}")
    
    return count_kept, count_discarded


def remove_adapters_with_cutadapt(
    r1_in: Path,
    r2_in: Path,
    r1_out: Path,
    r2_out: Path,
    threads: int = 24,
    logger: Optional[logging.Logger] = None
):
    """
    使用cutadapt去除接头序列
    
    参数:
        r1_in: 输入R1 FASTQ文件路径
        r2_in: 输入R2 FASTQ文件路径
        r1_out: 输出R1 FASTQ文件路径
        r2_out: 输出R2 FASTQ文件路径
        threads: 线程数（默认：24）
        logger: logger对象（可选）
    """
    if logger is None:
        logger = get_logger()
    
    logger.info(f"Removing adapters using cutadapt")
    logger.info(f"  Input R1: {r1_in}")
    logger.info(f"  Input R2: {r2_in}")
    logger.info(f"  Output R1: {r1_out}")
    logger.info(f"  Output R2: {r2_out}")
    logger.info(f"  Threads: {threads}")
    
    cutadapt = check_tool("cutadapt")
    
    # cutadapt命令：
    # -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC: R1 3'端接头
    # -A CCCCCCC: R2 3'端接头（7个C）
    # -m 20: 最小长度20bp
    # -j: 线程数
    cmd = (
        f"{cutadapt} "
        f"-j {threads} "
        f"-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC "
        f"-A CCCCCCC "
        f"-m 20 "
        f"-o {r1_out} "
        f"-p {r2_out} "
        f"{r1_in} "
        f"{r2_in}"
    )
    
    run_cmd(cmd, logger=logger, check=True)
    logger.info(f"Adapter removal completed: {r1_out}, {r2_out}")


def preprocess_paired_fastq(
    r1_in: Path,
    r2_in: Path,
    r1_out: Path,
    r2_out: Path,
    type: str,
    threads: int = 4,
    logger: Optional[logging.Logger] = None,
    tooldir: Optional[Path] = None
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
        tooldir: 工具目录（可选，用于日志路径）
    """
    if logger is None:
        # 如果没有提供logger，尝试创建logger
        if tooldir:
            log_dir = tooldir / "logs"
            log_dir.mkdir(parents=True, exist_ok=True)
            # 使用R1文件名按最后一个_分割后的[0]部分
            r1_prefix = r1_in.stem.rsplit('_', 1)[0] if '_' in r1_in.stem else r1_in.stem
            log_name = f"preprocess_{r1_prefix}.log"
            logger = setup_project_logger(log_dir, log_name)
        else:
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
        logger.info("RNA type: processing with new_library_preprocess-like workflow")
        # 固定参数：UMI 长度和 trim 长度、接头序列
        umi_length = 10
        trim_length = 3
        adapter_seq = "AAA"

        # 创建临时目录用于中间文件
        import tempfile
        temp_dir = r1_out.parent / ".tmp_preprocess_RNA"
        temp_dir.mkdir(parents=True, exist_ok=True)

        step1_r1 = temp_dir / "step1_r1_with_umi.fq"
        step2_r1 = temp_dir / "step2_dedup.fq"

        # 检查外部工具
        seqkit = check_tool("seqkit")
        cutadapt = check_tool("cutadapt")

        try:
            # Step 1: 按 new_library_preprocess 流程处理 R1/R2，生成带 UMI 的 R1（未去重）
            logger.info("RNA Step 1: 从 R2 选 11-13bp=AAA 的 pair，提取 R1/R2 UMI 并拼接到 R1 5' 端")
            total_pairs = 0
            kept_pairs = 0

            with open_read(r1_in) as r1_inf, open_read(r2_in) as r2_inf, open_write(step1_r1) as r1_out_step1:
                while True:
                    # 读取 R1 4 行
                    r1_header = r1_inf.readline()
                    if not r1_header:
                        break
                    r1_seq = r1_inf.readline().rstrip("\n")
                    r1_plus = r1_inf.readline()
                    r1_qual = r1_inf.readline().rstrip("\n")

                    # 读取 R2 4 行
                    r2_header = r2_inf.readline()
                    if not r2_header:
                        break
                    r2_seq = r2_inf.readline().rstrip("\n")
                    r2_plus = r2_inf.readline()
                    r2_qual = r2_inf.readline().rstrip("\n")

                    total_pairs += 1

                    # 检查 R2 的 11-13bp 是否为 AAA
                    if len(r2_seq) < 13 or r2_seq[10:13] != "AAA":
                        continue

                    # 提取 R2 和 R1 的前 10bp UMI
                    r2_umi = r2_seq[:umi_length]
                    r2_umi_qual = r2_qual[:umi_length]
                    r1_umi = r1_seq[:umi_length]
                    r1_umi_qual = r1_qual[:umi_length]

                    combined_umi = r1_umi + r2_umi
                    combined_umi_qual = r1_umi_qual + r2_umi_qual

                    # R1 去掉前 10bp 后的剩余部分
                    r1_remaining_seq = r1_seq[umi_length:]
                    r1_remaining_qual = r1_qual[umi_length:]

                    if not r1_remaining_seq:
                        continue

                    # 在 R1 剩余部分中，找到 adapter_seq 并截断（删除 adapter_seq 及之后所有内容）
                    if adapter_seq in r1_remaining_seq:
                        adapter_pos = r1_remaining_seq.find(adapter_seq)
                        trimmed_seq = r1_remaining_seq[:adapter_pos]
                        trimmed_qual = (
                            r1_remaining_qual[:adapter_pos]
                            if len(r1_remaining_qual) > adapter_pos
                            else r1_remaining_qual
                        )
                    else:
                        trimmed_seq = r1_remaining_seq
                        trimmed_qual = r1_remaining_qual

                    if not trimmed_seq:
                        continue

                    # 将 20bp UMI 拼接到 R1 5' 端
                    final_seq = combined_umi + trimmed_seq
                    final_qual = combined_umi_qual + trimmed_qual

                    # 质量值长度对齐
                    if len(final_qual) < len(final_seq):
                        final_qual = final_qual + "I" * (len(final_seq) - len(final_qual))
                    elif len(final_qual) > len(final_seq):
                        final_qual = final_qual[:len(final_seq)]

                    # 写入 Step1 R1
                    r1_out_step1.write(r1_header)
                    r1_out_step1.write(final_seq + "\n")
                    r1_out_step1.write(r1_plus)
                    r1_out_step1.write(final_qual + "\n")

                    kept_pairs += 1

            logger.info(f"RNA Step 1 completed: total_pairs={total_pairs}, kept_pairs={kept_pairs}")

            # Step 2: 使用 seqkit rmdup 去重
            logger.info("RNA Step 2: 使用 seqkit rmdup -s 去重")
            rmdup_cmd = f"{seqkit} rmdup -s -j {threads} {step1_r1} -o {step2_r1}"
            run_cmd(rmdup_cmd, logger=logger)

            # Step 3: 使用 cutadapt 去掉 UMI 和额外 trim 长度，保留长度 > 20 的序列
            total_trim = umi_length * 2 + trim_length  # 20bp UMI + 3bp trim_length = 23bp
            logger.info(
                f"RNA Step 3: 使用 cutadapt -u {total_trim} -m 20 去掉前 {total_trim}bp "
                f"(UMI 长度 {umi_length*2} + 额外 trim 长度 {trim_length})"
            )
            cut_cmd = (
                f"{cutadapt} -u {total_trim} -m 20 -j {threads} "
                f"-o {r1_out} {step2_r1}"
            )
            run_cmd(cut_cmd, logger=logger)

            # R2 在该新建库体系流程中不参与后续分析，这里简单复制原始 R2 以保持接口完整
            logger.info("RNA: Copying original R2 to output (R2 not modified in this workflow)")
            with open_read(r2_in) as inf, open_write(r2_out) as outf:
                shutil.copyfileobj(inf, outf)

            logger.info("RNA preprocessing completed")

        finally:
            # 清理临时文件
            try:
                if step1_r1.exists():
                    step1_r1.unlink()
                if step2_r1.exists():
                    step2_r1.unlink()
                if temp_dir.exists() and not any(temp_dir.iterdir()):
                    temp_dir.rmdir()
            except Exception as e:
                logger.warning(f"Failed to clean up RNA temporary files: {e}")

    elif type == 'DNA':
        logger.info("DNA type: processing R1 and R2")
        
        # 创建临时文件用于存储trimmed后的pair
        import tempfile
        temp_dir = r1_out.parent / ".tmp_preprocess"
        temp_dir.mkdir(parents=True, exist_ok=True)
        
        r1_trimmed = temp_dir / f"{r1_out.stem}.polyGtrim{r1_out.suffix}"
        r2_trimmed = temp_dir / f"{r2_out.stem}.polyGtrim{r2_out.suffix}"
        
        try:
            # Step 1: 切掉R1 5'端的所有G，根据G数量过滤pair
            logger.info("Step 1: Trimming 5' poly-G from R1 and filtering pairs")
            count_kept, count_discarded = trim_polyg_and_filter_pairs(
                r1_in=r1_in,
                r2_in=r2_in,
                r1_out=r1_trimmed,
                r2_out=r2_trimmed,
                min_g_count=7,
                logger=logger
            )
            
            # Step 2: 使用cutadapt去除接头
            logger.info("Step 2: Removing adapters using cutadapt")
            remove_adapters_with_cutadapt(
                r1_in=r1_trimmed,
                r2_in=r2_trimmed,
                r1_out=r1_out,
                r2_out=r2_out,
                threads=threads,
                logger=logger
            )
            
            logger.info("DNA preprocessing completed")
        finally:
            # 清理临时文件
            try:
                if r1_trimmed.exists():
                    r1_trimmed.unlink()
                if r2_trimmed.exists():
                    r2_trimmed.unlink()
                if temp_dir.exists() and not any(temp_dir.iterdir()):
                    temp_dir.rmdir()
            except Exception as e:
                logger.warning(f"Failed to clean up temporary files: {e}")
    
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
                       help="Log directory (optional, default: tooldir/logs)")
    parser.add_argument("--tooldir", type=Path, default=None,
                       help="Tool directory (for default log location)")
    
    args = parser.parse_args()
    
    # 设置日志
    if args.log_dir:
        log_dir = args.log_dir
    else:
        # 如果没有提供log_dir，使用tooldir/logs
        if args.tooldir:
            log_dir = args.tooldir / "logs"
        else:
            # 如果也没有tooldir，使用当前工作目录下的logs
            log_dir = Path.cwd() / "logs"
    
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # 生成日志文件名：preprocess_{r1按最后一个_分割后的[0]}.log
    # 例如：sample_001_R1.fastq.gz -> sample
    # setup_project_logger会自动在log_name后添加时间戳，最终格式为：
    # preprocess_{prefix}.log.{timestamp}
    r1_prefix = args.r1.stem.rsplit('_', 1)[0] if '_' in args.r1.stem else args.r1.stem
    log_name = f"preprocess_{r1_prefix}.log"
    
    logger = setup_project_logger(log_dir, log_name)
    
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
