#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
transTogenome.py

将转录组坐标的双端BAM文件转换为基因组坐标的BAM文件。

计算逻辑详解：
1. 转录组坐标到基因组坐标的转换：
   - 转录本坐标是连续的（0-based），不考虑内含子
   - 基因组坐标需要考虑外显子边界和内含子
   - 例如：转录本位置100可能对应基因组chr1:5000（如果前面有100bp的外显子）

2. CIGAR转换：
   - 转录组CIGAR：连续的M/D/I操作
   - 基因组CIGAR：需要在跨外显子边界时插入N（内含子）操作
   - 例如：转录组 "50M" 如果跨越两个外显子，可能变成 "30M 1000N 20M"

3. 链向处理：
   - 如果转录本是负链（strand='-'），需要：
     a. 对序列进行反向互补
     b. 对质量值进行反向
     c. 设置flag为16（反向链）
     d. CIGAR需要反向处理

4. 双端数据处理：
   - R1和R2需要分别转换
   - 需要正确设置mate信息（next_reference_id, next_reference_start等）
   - 需要计算template_length（插入片段长度）
   - 保持pair的一致性
   - 未配对的reads也会写入，但mate信息设置为未比对
"""

import pysam
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,FeatureLocation,CompoundLocation
from pathlib import Path
import logging
import re
import pandas as pd
from collections import defaultdict



def translate_cigar(
    old_cigar: list[tuple[int, int]],
    new_start: int,
    exons: list[tuple[int, int]],
    strand: str
) -> list[tuple[int, int]]:
    
    '''change cigar from transcriptome to genome'''
    
    if strand == "-":
        old_cigar = old_cigar[::-1]

    ref_pos = new_start
    exon_idx = 0
    exon_start, exon_end = exons[exon_idx]
    new_cigar = []

    for op, length in old_cigar:
        if op == 0 or op == 2:  # M or D
            remain = length
            while remain > 0:
                # If outside exon, jump to next exon and insert N
                while ref_pos >= exon_end:
                    exon_idx += 1
                    if exon_idx >= len(exons):
                        raise ValueError("Read exceeds exon boundaries")
                    next_start, next_end = exons[exon_idx]
                    intron_len = next_start - exon_end
                    if intron_len > 0:
                        new_cigar.append((3, intron_len))  # N
                    ref_pos = next_start
                    exon_start, exon_end = next_start, next_end

                chunk = min(remain, exon_end - ref_pos)
                if chunk > 0:
                    new_cigar.append((op, chunk))
                    ref_pos += chunk
                    remain -= chunk
                else:
                    # Prevent infinite loop
                    ref_pos = exon_end

        elif op == 1:  # I
            new_cigar.append((op, length))

        # Ignore H/S and other ops
        else:
            continue

    # Merge adjacent same ops, merge D+N or N+D to N
    merged = []
    i = 0
    while i < len(new_cigar):
        op, length = new_cigar[i]
        j = i + 1
        while j < len(new_cigar):
            next_op, next_len = new_cigar[j]
            if op == next_op:
                length += next_len
                j += 1
            elif (op == 2 and next_op == 3) or (op == 3 and next_op == 2):
                op = 3
                length += next_len
                j += 1
            else:
                break
        merged.append((op, length))
        i = j

    # Remove meaningless N at start/end (for safety)
    if merged and merged[0][0] == 3:
        merged = merged[1:]
    if merged and merged[-1][0] == 3:
        merged = merged[:-1]

    return merged



def change_segment(
    segment: pysam.AlignedSegment, 
    trans_dic: dict[str, dict],
    bam_out: pysam.AlignmentFile
) -> pysam.AlignedSegment:
    """
    将单个read从转录组坐标转换为基因组坐标
    
    计算逻辑：
    1. 提取转录本ID（从reference_name中，去掉后缀如_AGconvert）
    2. 查找转录本信息（染色体、链向、外显子坐标）
    3. 计算基因组起始位置：
       - 正链：从转录本起始位置累加外显子长度，找到对应的基因组位置
       - 负链：从转录本结束位置累加外显子长度，找到对应的基因组位置
    4. 转换CIGAR：在跨外显子边界时插入N（内含子）操作
    5. 处理链向：
       - 负链：序列反向互补，质量值反向，设置flag=16
    6. 处理双端数据的flag（paired, read1/read2, proper_pair等）
    注意：mate信息（next_reference_id, next_reference_start, template_length）将在主函数中设置
    
    参数:
        segment: 输入的AlignedSegment（转录组坐标）
        trans_dic: 转录本信息字典 {trans_id: {'chr': ..., 'strand': ..., 'starts': [...], 'ends': [...]}}
        bam_out: 输出BAM文件对象
    
    返回:
        转换后的AlignedSegment（基因组坐标），如果转换失败返回None
    """
    # 从reference_name中提取转录本ID：按|符号切割，取[0]
    old_chrname = segment.reference_name.split('|')[0]
    old_start = segment.reference_start
    old_end = segment.reference_end
    old_cigar = segment.cigartuples

    valid_ops = {0, 1, 2, 3}
    if not old_cigar or not all(op in valid_ops for op, _ in old_cigar):
        return None
    
    if old_chrname not in trans_dic:
        return None
    
    trans_info = trans_dic[old_chrname]
    new_chrname = trans_info['chr']
    current_length = 0
    new_strand = trans_info['strand']
    found = False
    
    # 计算基因组起始位置
    # 正链：从转录本起始位置（old_start）累加外显子长度
    # 负链：从转录本结束位置（old_end）累加外显子长度
    if new_strand == '+':
        for start, end in zip(trans_info['starts'], trans_info['ends']):
            exon_length = end - start + 1
            if current_length + exon_length <= old_start:
                current_length += exon_length
            else:
                # 找到对应的外显子，计算基因组位置
                new_start = start + old_start - current_length
                found = True
                break
    else:
        # 负链：从转录本结束位置计算
        for start, end in zip(trans_info['starts'], trans_info['ends']):
            exon_length = end - start + 1
            if current_length + exon_length <= old_end:
                current_length += exon_length
            else:
                # 负链时，old_end对应基因组位置
                new_start = start + old_end - current_length
                found = True
                break
    
    if not found:
        return None
    
    # 转换CIGAR
    exons = [(start, end) for start, end in zip(trans_info['starts'], trans_info['ends'])]
    try:
        new_cigar_tuple = translate_cigar(
            old_cigar=old_cigar, 
            new_start=new_start, 
            exons=exons,
            strand=new_strand
        )
    except ValueError:
        return None
    
    # 创建新的AlignedSegment
    segment_output = pysam.AlignedSegment()
    segment_output.tags = segment.tags.copy()

    qual = segment.query_qualities
    mpq = segment.mapping_quality
    seq = segment.query_sequence

    # 处理链向：负链需要反向互补序列和质量值
    if new_strand == "-":
        qual = qual[::-1] if qual is not None else None
        seq = str(Seq(seq).reverse_complement())
    
    # 注意：flag将在主函数中根据转录本链向和R1/R2正确设置
    # 这里先设置基本的flag，后续在主函数中会更新
    base_flag = 0
    if segment.is_paired:
        base_flag |= 0x1  # paired
        if segment.is_proper_pair:
            base_flag |= 0x2  # proper pair
        if segment.is_unmapped:
            base_flag |= 0x4  # unmapped
        if segment.mate_is_unmapped:
            base_flag |= 0x8  # mate unmapped
    
    segment_output.flag = base_flag
    segment_output.query_name = segment.query_name

    # 获取染色体ID
    ref_id = bam_out.get_tid(new_chrname)
    if ref_id is None:
        logging.warning(f"Chromosome {new_chrname} not found in BAM header. Available chromosomes: {[h['SN'] for h in bam_out.header['SQ']]}")
        return None
    segment_output.reference_id = ref_id

    segment_output.reference_start = new_start
    segment_output.cigartuples = new_cigar_tuple
    segment_output.query_sequence = seq
    segment_output.query_qualities = qual
    segment_output.mapping_quality = mpq
    
    # 注意：mate信息（next_reference_id, next_reference_start, template_length）
    # 将在主函数中处理，因为需要R1和R2都转换完成后才能正确设置
    
    return segment_output


def change_trans_to_genome(
    trans_bam: Path, 
    target_genome_fa: Path, 
    longest_transfile: Path, 
    outputbam: Path
) -> Path:
    """
    将转录组坐标的双端BAM文件转换为基因组坐标的BAM文件。
    
    计算流程：
    1. 读取转录本信息文件，构建转录本ID到基因组坐标的映射
    2. 读取输入BAM文件（转录组坐标，双端数据）
    3. 按read name分组，收集R1和R2
    4. 分别转换R1和R2到基因组坐标
    5. 计算mate信息（next_reference_id, next_reference_start, template_length）
    6. 只有R1和R2都转换成功时才写入输出BAM
    7. 未配对的reads也会写入，但mate信息设置为未比对
    
    参数:
        trans_bam: 输入BAM文件路径（转录组坐标，双端数据）
        target_genome_fa: 目标基因组FASTA文件路径
        longest_transfile: 最长转录本信息文件（TSV格式，包含chr, start, end, exon_starts, exon_ends, transcript_id等）
        outputbam: 输出BAM文件路径（基因组坐标）
    
    返回:
        输出BAM文件路径
    """
    bam_in = pysam.AlignmentFile(filename=trans_bam, mode='rb')

    logging.info('Reading transcriptome bam...')
    
    # 构建输出BAM的header
    header = {}
    header['HD'] = {'SO': 'unsorted', 'VN': '1.0'}
    header['SQ'] = []
    for seq in SeqIO.parse(target_genome_fa, 'fasta'):
        length = len(seq.seq)
        # 直接使用染色体名称，不添加后缀
        header['SQ'].append({"SN": seq.id, "LN": length})

    bam_out = pysam.AlignmentFile(outputbam, "wb", header=header)

    # 读取转录本信息
    # 格式：Chr\tmin_sites\tmax_sites\texon_starts\texon_ends\tTrans_length\tGene\tStrand\tTranscript_id
    trans_dic = {}  # pattern: {trans_id: {'chr': chr1, 'strand': '+', 'starts': [10,20,30], 'ends': [15,25,35]}}
    with open(longest_transfile, 'r') as file:
        firstline = file.readline()  # 跳过header
        for line in file:
            info = line.strip().split('\t')
            if len(info) < 9:
                continue
            transid = info[-1]  # Transcript_id在最后一列
            trans_dic[transid] = {
                'chr': info[0],  # 直接使用染色体名称，不添加后缀
                'strand': info[-2],  # Strand在倒数第二列
                'starts': [int(s) - 1 for s in info[3].split('_')],  # exon_starts，转为0-based
                'ends': [int(e) for e in info[4].split('_')]  # exon_ends，保持1-based
            }

    logging.info('Changing transcriptome coordinate to genome (paired-end only)')
    
    # 双端数据处理：按read name分组
    logging.info('Processing paired-end data...')
    pair_buffer = {}  # {read_name: {'r1': segment, 'r2': segment}}
    
    for segment in bam_in:
        read_name = segment.query_name
        
        # 初始化或更新pair信息
        if read_name not in pair_buffer:
            pair_buffer[read_name] = {}
        
        if segment.is_read1:
            pair_buffer[read_name]['r1'] = segment
        elif segment.is_read2:
            pair_buffer[read_name]['r2'] = segment
        else:
            # 如果不是R1或R2，跳过
            continue
        
        # 检查是否收集到完整的pair
        if 'r1' in pair_buffer[read_name] and 'r2' in pair_buffer[read_name]:
            r1 = pair_buffer[read_name]['r1']
            r2 = pair_buffer[read_name]['r2']
            
            # 转换R1
            r1_new = change_segment(
                segment=r1, 
                trans_dic=trans_dic,
                bam_out=bam_out
            )
            
            # 转换R2
            r2_new = change_segment(
                segment=r2, 
                trans_dic=trans_dic,
                bam_out=bam_out
            )
            
            # 只有R1和R2都转换成功时才写入
            if r1_new and r2_new:
                # 获取转录本链向信息（从R1或R2的reference_name中提取转录本ID）
                r1_trans_id = r1.reference_name.split('|')[0]
                if r1_trans_id in trans_dic:
                    transcript_strand = trans_dic[r1_trans_id]['strand']
                else:
                    # 如果找不到，使用R2的转录本ID
                    r2_trans_id = r2.reference_name.split('|')[0]
                    transcript_strand = trans_dic.get(r2_trans_id, {}).get('strand', '+')
                
                # 根据转录本链向设置flag
                # 正链转录本：R1在正链，R2在负链
                # 负链转录本：R1在负链，R2在正链
                if transcript_strand == '+':
                    # 正链转录本
                    # R1: 0x1 (paired) + 0x2 (proper) + 0x20 (mate reverse) + 0x40 (first) = 99
                    r1_new.flag = 0x1 | 0x2 | 0x20 | 0x40
                    # R2: 0x1 (paired) + 0x2 (proper) + 0x10 (reverse) + 0x80 (second) = 147
                    r2_new.flag = 0x1 | 0x2 | 0x10 | 0x80
                else:
                    # 负链转录本
                    # R1: 0x1 (paired) + 0x2 (proper) + 0x10 (reverse) + 0x40 (first) = 83
                    r1_new.flag = 0x1 | 0x2 | 0x10 | 0x40
                    # R2: 0x1 (paired) + 0x2 (proper) + 0x20 (mate reverse) + 0x80 (second) = 163
                    r2_new.flag = 0x1 | 0x2 | 0x20 | 0x80
                
                # 计算R1和R2的起始和结束位置（用于template_length计算）
                r1_start = r1_new.reference_start
                r1_end = r1_new.reference_start
                for op, length in r1_new.cigartuples:
                    if op == 0 or op == 2 or op == 3:  # M, D, N
                        r1_end += length
                
                r2_start = r2_new.reference_start
                r2_end = r2_new.reference_start
                for op, length in r2_new.cigartuples:
                    if op == 0 or op == 2 or op == 3:  # M, D, N
                        r2_end += length
                
                # 设置R1的mate信息
                r1_new.next_reference_id = r2_new.reference_id
                r1_new.next_reference_start = r2_new.reference_start
                
                # 设置R2的mate信息
                r2_new.next_reference_id = r1_new.reference_id
                r2_new.next_reference_start = r1_new.reference_start
                
                # 计算template_length（TLEN）
                if r1_new.reference_id == r2_new.reference_id and r1_new.reference_id >= 0:
                    # 同一染色体，计算template_length
                    left = min(r1_start, r2_start)
                    right = max(r1_end, r2_end)
                    tlen = right - left
                    
                    if r1_start == left:
                        r1_new.template_length = tlen
                        r2_new.template_length = -tlen
                    else:
                        r1_new.template_length = -tlen
                        r2_new.template_length = tlen
                else:
                    # 不同染色体或未比对
                    r1_new.template_length = 0
                    r2_new.template_length = 0
                
                # 写入R1和R2
                bam_out.write(r1_new)
                bam_out.write(r2_new)
            
            # 清除已处理的pair
            del pair_buffer[read_name]
    
    # 处理剩余的未配对reads（只写入R1或R2，但不设置mate信息）
    for read_name, pair_data in pair_buffer.items():
        if 'r1' in pair_data:
            r1_new = change_segment(
                segment=pair_data['r1'], 
                trans_dic=trans_dic,
                bam_out=bam_out
            )
            if r1_new:
                # 获取转录本链向
                r1_trans_id = pair_data['r1'].reference_name.split('|')[0]
                transcript_strand = trans_dic.get(r1_trans_id, {}).get('strand', '+')
                
                # 设置未配对R1的flag
                if transcript_strand == '+':
                    # 正链转录本，R1在正链
                    r1_new.flag = 0x1 | 0x8 | 0x40  # paired + mate unmapped + first
                else:
                    # 负链转录本，R1在负链
                    r1_new.flag = 0x1 | 0x8 | 0x10 | 0x40  # paired + mate unmapped + reverse + first
                
                # 未配对的read，设置mate信息为未比对
                r1_new.next_reference_id = -1
                r1_new.next_reference_start = -1
                r1_new.template_length = 0
                bam_out.write(r1_new)
        elif 'r2' in pair_data:
            r2_new = change_segment(
                segment=pair_data['r2'], 
                trans_dic=trans_dic,
                bam_out=bam_out
            )
            if r2_new:
                # 获取转录本链向
                r2_trans_id = pair_data['r2'].reference_name.split('|')[0]
                transcript_strand = trans_dic.get(r2_trans_id, {}).get('strand', '+')
                
                # 设置未配对R2的flag
                if transcript_strand == '+':
                    # 正链转录本，R2在负链
                    r2_new.flag = 0x1 | 0x8 | 0x10 | 0x80  # paired + mate unmapped + reverse + second
                else:
                    # 负链转录本，R2在正链
                    r2_new.flag = 0x1 | 0x8 | 0x80  # paired + mate unmapped + second
                
                # 未配对的read，设置mate信息为未比对
                r2_new.next_reference_id = -1
                r2_new.next_reference_start = -1
                r2_new.template_length = 0
                bam_out.write(r2_new)

    bam_in.close()
    bam_out.close()
    logging.info('Finished changing!')
    
    return outputbam

if __name__ == "__main__":
    trans_bam = Path(r"test/test.target.bam")
    target_genome_fa = Path(r"raw_data/human/GRCh38.p14.genome.fa")
    longest_transfile = Path(r"raw_data/human/gencode.v49.transcripts.longest.trans")
    outputbam = Path(r"test/test.target_genome.bam")
    change_trans_to_genome(trans_bam, target_genome_fa, longest_transfile, outputbam)