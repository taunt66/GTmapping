#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
filter_and_return_origin_bam.py

根据R1和R2的碱基索引文件将BAM文件恢复到原始状态，并过滤含有过量指定碱基的pair。

功能：
1. 根据索引文件将指定位置的碱基替换回原始碱基
2. 过滤功能：
   - 指定碱基（如A）和cutoff（如3）
   - 对于正链reads（无论R1还是R2）：统计指定碱基的个数
   - 对于负链reads（无论R1还是R2）：统计互补碱基的个数（A->T, C->G等）
   - 每个pair必须R1和R2都小于cutoff才保留
3. 输出两个BAM文件：
   - origin_bam：恢复到原始的BAM
   - filter_bam：过滤后的BAM
"""

import pysam
from pathlib import Path
from typing import Optional, Tuple
import logging
import re
from collections import defaultdict

# 互补碱基映射表
REVERSE_COMPLEMENT_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}


def parse_index_file(index_file: Path) -> dict:
    """
    解析index文件，返回 {read_id: {base: [indices]}} 字典
    
    例如：{'read1': {'A': [36, 66], 'C': [17]}, 'read2': {'A': [], 'C': []}}
    
    参数:
        index_file: index文件路径
    
    返回:
        dict: {read_id: {base: [indices]}}
    """
    index_dict = {}
    
    with open(index_file, 'r') as f:
        # 读取标题行
        header = f.readline().strip()
        header_fields = header.split('\t')
        
        # 第一列是read_id，从第二列开始查找所有*_index列
        base_to_col = {}  # {'A': 1, 'C': 2, ...}
        
        for i, col in enumerate(header_fields[1:], start=1):  # 从第二列开始
            col_upper = col.upper()
            # 匹配以_base_index结尾的列名
            match = re.match(r'^([ACGT])_INDEX$', col_upper)
            if match:
                base = match.group(1)
                base_to_col[base] = i
        
        if not base_to_col:
            raise ValueError(f"在index文件中未找到任何*_index列（如A_index, C_index等）")
        
        # 读取数据
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue
            
            read_id = fields[0]
            read_indices = {}
            
            # 解析每个碱基类型的索引
            for base, col_idx in base_to_col.items():
                if len(fields) > col_idx:
                    indices_str = fields[col_idx]
                    if indices_str.upper() == 'NA' or indices_str == '':
                        read_indices[base] = []
                    else:
                        # 解析索引（支持用_分隔的多个索引）
                        try:
                            indices = [int(x) for x in indices_str.split('_')]
                            read_indices[base] = indices
                        except ValueError:
                            logging.warning(f"无法解析索引 '{indices_str}' for read {read_id}, 跳过")
                            read_indices[base] = []
                else:
                    read_indices[base] = []
            
            index_dict[read_id] = read_indices
    
    return index_dict


def restore_sequence(read: pysam.AlignedSegment, read_indices: dict) -> bool:
    """
    根据索引文件恢复read的序列
    
    参数:
        read: pysam AlignedSegment对象
        read_indices: {base: [indices]} 字典
    
    返回:
        bool: 是否修改了序列
    """
    seq = list(read.query_sequence)
    read_len = len(seq)
    is_reverse = read.is_reverse
    modified = False
    
    # 遍历所有碱基类型（A, C等）
    for base, indices in read_indices.items():
        if indices:  # 如果有索引需要替换
            for i in indices:
                if 0 <= i < read_len:
                    if is_reverse:
                        # 对于反向的read，需要：
                        # 1. 坐标转换：原始坐标i在反向互补序列中的位置是 read_len - 1 - i
                        # 2. 碱基互补：原始碱基base需要转换为互补碱基
                        ri = read_len - 1 - i
                        if 0 <= ri < read_len:
                            seq[ri] = REVERSE_COMPLEMENT_DICT[base]
                            modified = True
                    else:
                        # 正向read，直接替换
                        seq[i] = base
                        modified = True
    
    if modified:
        read.query_sequence = ''.join(seq)
    
    return modified


def count_base_in_read(read: pysam.AlignedSegment, target_base: str) -> int:
    """
    统计read中指定碱基的个数（考虑链向）
    
    参数:
        read: pysam AlignedSegment对象
        target_base: 要统计的碱基（如'A'）
    
    返回:
        int: 碱基个数
    """
    seq = read.query_sequence or ""
    
    if read.is_reverse:
        # 负链：统计互补碱基
        complement_base = REVERSE_COMPLEMENT_DICT[target_base]
        count = seq.count(complement_base) + seq.count(complement_base.lower())
    else:
        # 正链：统计原始碱基
        count = seq.count(target_base) + seq.count(target_base.lower())
    
    return count


def filter_and_return_origin_bam(
    bam_path: Path,
    r1_index_file: Path,
    r2_index_file: Path,
    origin_bam: Path,
    filter_bam: Path,
    target_base: str = 'A',
    cutoff: int = 3,
    logger: Optional[logging.Logger] = None
) -> Tuple[Path, Path]:
    """
    根据索引文件恢复BAM到原始状态，并过滤含有过量指定碱基的pair
    
    参数:
        bam_path: 输入BAM文件路径（双端数据）
        r1_index_file: R1的索引文件路径
        r2_index_file: R2的索引文件路径
        origin_bam: 输出原始BAM文件路径（恢复后的BAM）
        filter_bam: 输出过滤后的BAM文件路径
        target_base: 要检查的碱基（默认'A'）
        cutoff: 碱基个数阈值（默认3），R1和R2都必须小于此值才保留
        logger: logger对象（可选）
    
    返回:
        (origin_bam_path, filter_bam_path): 输出的两个BAM文件路径的元组
    """
    if logger is None:
        logger = logging.getLogger(__name__)
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('[%(levelname)s] %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
    
    # 验证target_base
    if target_base.upper() not in ['A', 'C', 'G', 'T']:
        raise ValueError(f"target_base必须是A、C、G或T之一，当前值: {target_base}")
    target_base = target_base.upper()
    
    # 确保输出目录存在
    origin_bam.parent.mkdir(parents=True, exist_ok=True)
    filter_bam.parent.mkdir(parents=True, exist_ok=True)
    
    # 读取索引文件
    logger.info(f"读取R1索引文件: {r1_index_file}")
    r1_index_dict = parse_index_file(r1_index_file)
    logger.info(f"R1索引记录数: {len(r1_index_dict)}")
    
    logger.info(f"读取R2索引文件: {r2_index_file}")
    r2_index_dict = parse_index_file(r2_index_file)
    logger.info(f"R2索引记录数: {len(r2_index_dict)}")
    
    # 打开输入和输出BAM文件
    bam_in = pysam.AlignmentFile(bam_path, 'rb')
    bam_origin = pysam.AlignmentFile(origin_bam, 'wb', template=bam_in)
    bam_filter = pysam.AlignmentFile(filter_bam, 'wb', template=bam_in)
    
    # 按read name分组处理pair
    # 注意：pysam的迭代器返回的是同一个对象，需要保存序列的副本
    pair_buffer = {}  # {read_name: {'r1': read, 'r2': read}}
    
    processed_pairs = 0
    restored_pairs = 0
    filtered_pairs = 0
    kept_pairs = 0
    
    for read in bam_in:
        read_name = read.query_name
        
        # 初始化或更新pair信息
        if read_name not in pair_buffer:
            pair_buffer[read_name] = {}
        
        # 创建read的副本（保存序列和其他属性）
        read_copy = pysam.AlignedSegment()
        read_copy.query_name = read.query_name
        read_copy.flag = read.flag
        read_copy.reference_id = read.reference_id
        read_copy.reference_start = read.reference_start
        read_copy.mapping_quality = read.mapping_quality
        read_copy.cigartuples = read.cigartuples
        read_copy.next_reference_id = read.next_reference_id
        read_copy.next_reference_start = read.next_reference_start
        read_copy.template_length = read.template_length
        read_copy.query_sequence = read.query_sequence
        read_copy.query_qualities = read.query_qualities
        read_copy.tags = read.tags
        
        if read.is_read1:
            pair_buffer[read_name]['r1'] = read_copy
        elif read.is_read2:
            pair_buffer[read_name]['r2'] = read_copy
        else:
            # 如果不是R1或R2，跳过
            continue
        
        # 检查是否收集到完整的pair
        if 'r1' in pair_buffer[read_name] and 'r2' in pair_buffer[read_name]:
            r1 = pair_buffer[read_name]['r1']
            r2 = pair_buffer[read_name]['r2']
            
            processed_pairs += 1
            
            # 恢复R1序列
            r1_modified = False
            if read_name in r1_index_dict:
                r1_modified = restore_sequence(r1, r1_index_dict[read_name])
            
            # 恢复R2序列
            r2_modified = False
            if read_name in r2_index_dict:
                r2_modified = restore_sequence(r2, r2_index_dict[read_name])
            
            if r1_modified or r2_modified:
                restored_pairs += 1
            
            # 写入原始BAM（恢复后的）
            bam_origin.write(r1)
            bam_origin.write(r2)
            
            # 统计碱基个数并过滤
            r1_count = count_base_in_read(r1, target_base)
            r2_count = count_base_in_read(r2, target_base)
            
            # 如果R1和R2都小于cutoff，保留这个pair
            if r1_count < cutoff and r2_count < cutoff:
                bam_filter.write(r1)
                bam_filter.write(r2)
                kept_pairs += 1
            else:
                filtered_pairs += 1
            
            # 清除已处理的pair
            del pair_buffer[read_name]
    
    # 处理剩余的未配对reads
    for read_name, pair_data in pair_buffer.items():
        if 'r1' in pair_data:
            r1 = pair_data['r1']
            if read_name in r1_index_dict:
                restore_sequence(r1, r1_index_dict[read_name])
            bam_origin.write(r1)
            # 未配对的read也进行过滤
            r1_count = count_base_in_read(r1, target_base)
            if r1_count < cutoff:
                bam_filter.write(r1)
        elif 'r2' in pair_data:
            r2 = pair_data['r2']
            if read_name in r2_index_dict:
                restore_sequence(r2, r2_index_dict[read_name])
            bam_origin.write(r2)
            # 未配对的read也进行过滤
            r2_count = count_base_in_read(r2, target_base)
            if r2_count < cutoff:
                bam_filter.write(r2)
    
    bam_in.close()
    bam_origin.close()
    bam_filter.close()
    
    logger.info(f"处理完成:")
    logger.info(f"  处理的pair数: {processed_pairs}")
    logger.info(f"  恢复的pair数: {restored_pairs}")
    logger.info(f"  保留的pair数: {kept_pairs}")
    logger.info(f"  过滤的pair数: {filtered_pairs}")
    logger.info(f"  原始BAM: {origin_bam}")
    logger.info(f"  过滤BAM: {filter_bam}")
    
    return origin_bam, filter_bam

