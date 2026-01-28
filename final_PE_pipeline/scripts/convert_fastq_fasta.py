#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
convert_fastq_fasta.py

提供两个函数用于碱基转换：
1. change_fastq: 转换FASTQ文件，支持gz压缩和纯文本格式
2. change_fasta: 转换FASTA文件

两个函数都支持多个转换规则（如 A:G, C:T），返回值都是生成文件的Path对象。
"""

from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
import re
from typing import List, Tuple, Optional


def _open_read_handle(path: Path):
    """根据文件扩展名自动选择读取方式（支持gz压缩）"""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    else:
        return open(path, "r")


def _open_write_handle(path: Path):
    """根据文件扩展名自动选择写入方式（支持gz压缩）"""
    if str(path).endswith(".gz"):
        # 使用compresslevel=6作为默认压缩级别（平衡速度和压缩率）
        return gzip.open(path, "wt", compresslevel=6)
    else:
        return open(path, "w")


def _parse_convert_pairs(convert_list: List[str]) -> Tuple[dict, List[str], str]:
    """
    解析转换规则列表
    
    参数:
        convert_list: 转换规则列表，如 ['A:G', 'C:T']
    
    返回:
        (mapping_dict, pres_list, label)
        - mapping_dict: 碱基映射字典，如 {'A': 'G', 'C': 'T'}
        - pres_list: 原始碱基列表，如 ['A', 'C']
        - label: 标签字符串，如 'A-G_C-T'
    """
    mapping = {}
    pres = []
    parts = []
    for item in convert_list:
        if ":" not in item:
            raise ValueError(f"Invalid convert spec '{item}'. Expect form X:Y")
        pre, post = item.split(":", 1)
        pre = pre.strip().upper()
        post = post.strip().upper()
        if len(pre) != 1 or len(post) != 1 or pre not in "ACGT" or post not in "ACGT":
            raise ValueError(f"Invalid bases in convert '{item}'. Use single A/C/G/T.")
        if pre in mapping:
            raise ValueError(f"Duplicate mapping for base {pre}")
        mapping[pre] = post
        pres.append(pre)
        parts.append(f"{pre}-{post}")
    label = "_".join(parts)
    return mapping, pres, label


def change_fastq(input_fastq: Path, output_dir: Path, 
                 convert_list: List[str], 
                 output_fastq: Optional[Path] = None) -> Tuple[Path, Path]:
    """
    转换FASTQ文件，支持gz压缩和纯文本格式
    
    参数:
        input_fastq: 输入FASTQ文件路径（支持.fq, .fastq, .fq.gz, .fastq.gz）
        output_dir: 输出目录
        convert_list: 转换规则列表，如 ['A:G', 'C:T']
        output_fastq: 可选，指定输出FASTQ文件路径（如果指定，将使用此路径）
    
    返回:
        (changed_fastq_path, index_file_path)
        - changed_fastq_path: 转换后的FASTQ文件路径
        - index_file_path: 原始碱基位置索引文件路径
    """
    # 解析转换规则
    mapping, pres, label = _parse_convert_pairs(convert_list)
    
    # 确保输出目录存在
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 如果指定了输出文件路径，直接使用
    if output_fastq is not None:
        changed_fq_file = output_fastq
        # 根据输出文件名生成索引文件名
        base_name = output_fastq.stem
        if base_name.endswith(".fq") or base_name.endswith(".fastq"):
            base_name = base_name.rsplit(".", 1)[0]
        index_file = output_dir / f"{base_name}_index_{label}.txt"
    else:
        # 确定输出文件名（保持输入格式）
        input_stem = input_fastq.stem
        # 如果输入是.gz，stem会去掉.gz，需要手动处理
        if str(input_fastq).endswith(".gz"):
            # 如果stem是.fastq或.fq，需要保留
            if input_stem.endswith(".fastq") or input_stem.endswith(".fq"):
                base_name = input_stem
            else:
                base_name = input_fastq.name.replace(".gz", "")
        else:
            base_name = input_fastq.name
        
        # 确定输出文件扩展名
        if base_name.endswith(".fastq"):
            ext = ".fastq"
            base_name = base_name.replace(".fastq", "")
        elif base_name.endswith(".fq"):
            ext = ".fq"
            base_name = base_name.replace(".fq", "")
        else:
            ext = ".fastq"
        
        # 构建输出文件名（根据输入格式决定输出格式）
        out_name = f"{base_name}_{label}_changed{ext}"
        input_is_gz = str(input_fastq).endswith(".gz")
        if input_is_gz:
            out_name += ".gz"
        
        changed_fq_file = output_dir / out_name
        index_file = output_dir / f"{base_name}_index_{label}.txt"
    
    # 处理文件
    with (_open_read_handle(input_fastq) as infh, 
          _open_write_handle(changed_fq_file) as outfh, 
          _open_write_handle(index_file) as idxfh):
        # 写入索引文件标题
        header = ["read_id"] + [f"{p}_index" for p in pres]
        idxfh.write("\t".join(header) + "\n")
        
        for record in SeqIO.parse(infh, "fastq"):
            seq = str(record.seq)
            seq_upper = seq.upper()
            
            # 计算每个原始碱基的位置索引（基于原始序列）
            idx_cols = []
            for pre in pres:
                positions = [str(m.start()) for m in re.finditer(pre, seq_upper)]
                if positions:
                    idx_cols.append("_".join(positions))
                else:
                    idx_cols.append("NA")
            
            # 基于原始序列同时应用所有转换
            out_chars = []
            for ch in seq_upper:
                if ch in mapping:
                    out_chars.append(mapping[ch])
                else:
                    out_chars.append(ch)
            new_seq = "".join(out_chars)
            
            # 保留质量分数
            qual = record.letter_annotations.get("phred_quality", None)
            
            # 创建新的SeqRecord
            new_rec = SeqRecord(
                Seq(new_seq),
                id=record.id,
                name=record.name,
                description=record.description
            )
            if qual is not None:
                new_rec.letter_annotations["phred_quality"] = qual
            
            SeqIO.write(new_rec, outfh, "fastq")
            idxfh.write(record.id + "\t" + "\t".join(idx_cols) + "\n")
    
    return changed_fq_file, index_file


def change_fasta(input_fasta: Path, output_fasta: Path, 
                 convert_list: List[str]) -> Path:
    """
    转换FASTA文件
    
    参数:
        input_fasta: 输入FASTA文件路径
        output_fasta: 输出FASTA文件路径
        convert_list: 转换规则列表，如 ['A:G', 'C:T']
    
    返回:
        输出FASTA文件路径（Path对象）
    """
    # 解析转换规则
    mapping, pres, label = _parse_convert_pairs(convert_list)
    
    # 确保输出目录存在
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    
    # 处理文件
    with open(input_fasta, 'r') as infh, open(output_fasta, 'w') as outfh:
        for record in SeqIO.parse(infh, format='fasta'):
            seq_upper = str(record.seq).upper()
            
            # 基于原始序列同时应用所有转换
            out_chars = []
            for ch in seq_upper:
                if ch in mapping:
                    out_chars.append(mapping[ch])
                else:
                    out_chars.append(ch)
            new_seq = Seq("".join(out_chars))
            
            # 保持原始ID不变
            new_record = SeqRecord(
                seq=new_seq,
                id=record.id,
                name=record.name,
                description=record.description
            )
            SeqIO.write(new_record, handle=outfh, format='fasta')
    
    return output_fasta


if __name__ == "__main__":
    input_fasta1=Path(r"/home/wangjing_pkuhpc/lustre2/GT_mapping_pipeline/raw_data/human/GRCh38.p14.genome.fa")
    input_fasta2=Path(r"/home/wangjing_pkuhpc/lustre2/GT_mapping_pipeline/raw_data/human/gencode.v49.transcripts.fa")
    output_fasta1=Path(r"/home/wangjing_pkuhpc/lustre2/GT_mapping_pipeline/raw_data/human/GRCh38.p14.genome.fa_A-G_C-T_changed.fasta")
    output_fasta2=Path(r"/home/wangjing_pkuhpc/lustre2/GT_mapping_pipeline/raw_data/human/gencode.v49.transcripts.fa_A-G_C-T_changed.fasta")
    output_fasta3=Path(r"/home/wangjing_pkuhpc/lustre2/GT_mapping_pipeline/raw_data/human/GRCh38.p14.genome.fa_T_C_G_A_changed.fasta")
    
    convert_list=['A:G', 'C:T']
    convert_list2=['T:C', 'G:A']
    change_fasta(input_fasta1, output_fasta1, convert_list)
    change_fasta(input_fasta2, output_fasta2, convert_list)
    change_fasta(input_fasta1, output_fasta3, convert_list2)