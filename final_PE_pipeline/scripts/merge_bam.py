#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
merge_bam.py

合并多个基因组坐标的双端BAM文件。
处理不同BAM文件中染色体名称的差异（如chr1_AGconvert -> chr1），统一为规范化的染色体名称。
"""

import pysam
from pathlib import Path
import shutil
import logging
from typing import List, Optional

# 导入日志和命令运行函数
try:
    from .run_cmd_and_logging import run_cmd, get_logger
except ImportError:
    from run_cmd_and_logging import run_cmd, get_logger


def read_headers(fn: str, hid: int, new_header: dict, hid_dict: dict, lift_over: dict):
    """
    读取BAM文件的header，构建染色体名称映射
    
    参数:
        fn: BAM文件路径
        hid: 当前header ID计数器
        new_header: 新的BAM header字典
        hid_dict: 规范化染色体名称到header ID的映射
        lift_over: 旧BAM的tid到新BAM tid的映射
    
    返回:
        (hid, new_header, hid_dict, lift_over)
    """
    lift_over[fn] = {}
    with pysam.AlignmentFile(fn, 'rb') as INPUT:
        n = 0
        for header in INPUT.header['SQ']:
            raw_name = header['SN']
            # 规范化染色体名称：去掉后缀（如chr1_AGconvert -> chr1）
            canonical_name = raw_name.rsplit('_', 1)[0]
            if canonical_name not in hid_dict:
                hid_dict[canonical_name] = hid
                new_header['SQ'].append({'SN': canonical_name, 'LN': header['LN']})
                hid += 1
            lift_over[fn][n] = hid_dict[canonical_name]  # 旧 bam 的 tid 映射为新 bam 中 canonical tid
            n += 1
    return hid, new_header, hid_dict, lift_over


def process_and_write_bam(bam_path: Path, output_bam: Path, lift_over: dict, new_header: dict):
    """
    处理单个BAM文件，映射reference_id并写入临时BAM文件
    
    参数:
        bam_path: 输入BAM文件路径
        output_bam: 输出BAM文件路径
        lift_over: reference_id映射字典
        new_header: 新的BAM header
    """
    bam_in = pysam.AlignmentFile(bam_path, 'rb')
    bam_out = pysam.AlignmentFile(output_bam, 'wb', header=new_header)

    for read in bam_in:
        # 替换 reference_id（映射到规范化的染色体ID）
        old_tid = read.reference_id
        try:
            read.reference_id = lift_over[str(bam_path)][old_tid]
        except KeyError:
            continue  # 忽略在 header 中不存在的 reference

        bam_out.write(read)

    bam_in.close()
    bam_out.close()


def merge_bams(
    bams: List[Path], 
    outputdir: Path, 
    prefix: str,
    logger: Optional[logging.Logger] = None
) -> Path:
    """
    合并多个基因组坐标的双端BAM文件
    
    参数:
        bams: 输入BAM文件路径列表
        outputdir: 输出目录
        prefix: 输出文件前缀
        logger: logger对象（可选）
    
    返回:
        合并并排序后的BAM文件路径
    """
    if logger is None:
        logger = get_logger()
    
    outputdir.mkdir(parents=True, exist_ok=True)

    temp_dir = outputdir / '.tmp_mergebam'
    temp_dir.mkdir(parents=True, exist_ok=True)

    new_header = {'HD': {'VN': '1.0', 'SO': 'coordinate'}, 'SQ': []}
    hid_dict = {}
    lift_over = {}
    hid = 0

    # 读取所有BAM文件的header，构建统一的header
    for bam in bams:
        logger.info(f'Reading header from {bam}')
        hid, new_header, hid_dict, lift_over = read_headers(str(bam), hid, new_header, hid_dict, lift_over)

    # 处理每个BAM文件，映射reference_id
    temp_bams = []
    for i, bam in enumerate(bams):
        logger.info(f'Processing BAM {i+1}/{len(bams)}: {bam}')
        temp_bam = temp_dir / f"temp_{i}.bam"
        process_and_write_bam(bam, temp_bam, lift_over, new_header)
        temp_bams.append(str(temp_bam))

    # 合并所有临时BAM文件
    merged_bam = outputdir / f"{prefix}_merged.bam"
    sorted_bam = outputdir / f"{prefix}_merged.sorted.bam"
    
    # 使用samtools merge（通过run_cmd）
    logger.info(f'Merging BAM files...')
    temp_bams_str = " ".join(temp_bams)
    merge_cmd = f"samtools merge -f {merged_bam} {temp_bams_str}"
    run_cmd(merge_cmd, logger=logger)
    
    # 使用samtools sort（通过run_cmd）
    logger.info(f'Sorting merged BAM...')
    sort_cmd = f"samtools sort -o {sorted_bam} {merged_bam}"
    run_cmd(sort_cmd, logger=logger)
    
    # 使用samtools index（通过run_cmd）
    logger.info(f'Indexing sorted BAM...')
    index_cmd = f"samtools index {sorted_bam}"
    run_cmd(index_cmd, logger=logger)

    # 清理临时文件
    shutil.rmtree(temp_dir)
    logger.info(f"BAM merging completed: {sorted_bam}")
    return sorted_bam


