#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_target_aln&dump_others.py

从双端 BAM 文件中提取符合条件的配对（R1 和 R2 都是 unique primary alignment，且 R1 满足链向要求），
并将未被选到的 pair（包括 unmapped、链向不满足要求、多重比对的）输出到 FASTQ 文件。
"""

import pysam
from pathlib import Path
from collections import defaultdict
from typing import Optional, Tuple
import gzip
import subprocess
import shutil
import logging

# 处理相对导入和绝对导入
try:
    from .run_cmd_and_logging import run_cmd, get_logger
except ImportError:
    # 如果相对导入失败，尝试绝对导入
    import sys
    from pathlib import Path as P
    
    # 添加scripts目录到路径
    script_dir = P(__file__).parent
    if str(script_dir) not in sys.path:
        sys.path.insert(0, str(script_dir))
    
    try:
        from run_cmd_and_logging import run_cmd, get_logger
    except ImportError:
        # 如果都失败，创建一个简单的fallback
        def get_logger():
            logger = logging.getLogger("pipeline")
            if not logger.handlers:
                handler = logging.StreamHandler()
                formatter = logging.Formatter('[%(levelname)s] %(message)s')
                handler.setFormatter(formatter)
                logger.addHandler(handler)
                logger.setLevel(logging.INFO)
            return logger
        
        def run_cmd(cmd, logger=None, **kwargs):
            if logger is None:
                logger = get_logger()
            logger.info(f"[CMD] {cmd}")
            result = subprocess.run(cmd, shell=True, **kwargs)
            if result.returncode != 0:
                logger.error(f"[ERROR] Command failed with return code {result.returncode}")
            return result


def is_primary_unique(aln: pysam.AlignedSegment) -> bool:
    """
    检查 alignment 是否是 primary 且 unique。
    
    判断逻辑：
    1. 优先根据 NH tag 判断：
       - 如果有 NH tag，则 NH == 1 表示 unique
    2. 如果没有 NH tag，则根据其他标志判断：
       - 如果是 secondary、supplementary 或 unmapped，则不是 unique
       - 否则认为是 unique

    Args:
        aln: pysam AlignedSegment

    Returns:
        bool: True 如果是 primary 且 unique
    """
    # 优先检查 NH tag
    if aln.has_tag("NH"):
        nh = aln.get_tag("NH")
        # NH == 1 表示 unique alignment
        return nh == 1
    
    # 如果没有 NH tag，根据其他标志判断
    # 如果是 secondary、supplementary 或 unmapped，则不是 unique
    if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
        return False
    
    # 否则认为是 unique（primary alignment 且没有多重比对标志）
    return True


def _open_fastq(file_path: Path, mode: str = 'r'):
    """打开 FASTQ 文件，自动处理 gzip 压缩"""
    if file_path.suffix == '.gz':
        return gzip.open(file_path, mode + 't')
    else:
        return open(file_path, mode)


def _normalize_read_name(header_line: str) -> str:
    """标准化 read name，去掉 @ 符号和 /1、/2 后缀"""
    name = header_line.strip().split()[0].lstrip('@')
    # 去掉 /1 或 /2 后缀
    if name.endswith('/1') or name.endswith('/2'):
        name = name[:-2]
    return name


def _read_fastq_pair(r1_file, r2_file):
    """读取一对 FASTQ 记录（R1 和 R2）"""
    r1_header = r1_file.readline()
    if not r1_header:
        return None, None
    
    r1_seq = r1_file.readline().rstrip('\n')
    r1_plus = r1_file.readline()
    r1_qual = r1_file.readline().rstrip('\n')
    
    r2_header = r2_file.readline()
    if not r2_header:
        return None, None
    
    r2_seq = r2_file.readline().rstrip('\n')
    r2_plus = r2_file.readline()
    r2_qual = r2_file.readline().rstrip('\n')
    
    # 标准化 read name
    r1_name = _normalize_read_name(r1_header)
    r2_name = _normalize_read_name(r2_header)
    
    r1_record = (r1_header, r1_seq, r1_plus, r1_qual, r1_name)
    r2_record = (r2_header, r2_seq, r2_plus, r2_qual, r2_name)
    
    return r1_record, r2_record


def extract_target_aln_and_dump_others(
    bam_path: Path,
    output_bam: Path,
    output_r1_fastq: Path,
    output_r2_fastq: Path,
    r1_strand: Optional[str] = None,  # 'plus', 'minus', or None (all)
    logger: Optional[logging.Logger] = None
) -> Tuple[Path, Path, Path]:
    """
    从双端 BAM 文件中提取符合条件的配对，并将未选中的 pair 输出到 FASTQ 文件。
    
    参数:
        bam_path: 输入 BAM 文件路径
        output_bam: 输出 BAM 文件路径（符合条件的 pair）
        output_r1_fastq: 输出 R1 FASTQ 文件路径（未选中的 pair，排序前）
        output_r2_fastq: 输出 R2 FASTQ 文件路径（未选中的 pair，排序前）
        r1_strand: R1 的链向要求，'plus'、'minus' 或 None（提取所有方向）
        logger: logger对象（可选）
    
    返回:
        (output_bam_path, repaired_r1_fastq_path, repaired_r2_fastq_path): 
        输出的 BAM 和排序后的 R1、R2 FASTQ 文件路径的元组
    """
    if logger is None:
        logger = get_logger()
    
    if r1_strand is not None and r1_strand not in ["plus", "minus"]:
        raise ValueError("r1_strand must be 'plus', 'minus', or None")
    
    # 确保输出目录存在
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    output_r1_fastq.parent.mkdir(parents=True, exist_ok=True)
    output_r2_fastq.parent.mkdir(parents=True, exist_ok=True)
    
    # 打开输入 BAM
    in_bam = pysam.AlignmentFile(str(bam_path), "rb")
    
    # 存储每个 read name 的 R1 和 R2 alignment
    # read_dict[read_name] = {'r1': [aln1, aln2, ...], 'r2': [aln1, aln2, ...]}
    read_dict = defaultdict(lambda: {"r1": [], "r2": []})
    
    # 第一遍扫描：收集所有 primary unique 的 R1 和 R2
    for aln in in_bam:
        # 只处理 primary unique alignment
        if not is_primary_unique(aln):
            continue
        
        read_name = aln.query_name
        
        if aln.is_read1:
            read_dict[read_name]["r1"].append(aln)
        elif aln.is_read2:
            read_dict[read_name]["r2"].append(aln)
    
    # 第二遍：筛选符合条件的配对，并保存 alignment 对象
    valid_pairs = {}  # 存储符合条件的 pair: {read_name: (r1_aln, r2_aln)}
    
    for read_name, reads in read_dict.items():
        r1_list = reads["r1"]
        r2_list = reads["r2"]
        
        # 必须同时有 R1 和 R2
        if not r1_list or not r2_list:
            continue
        
        # 如果 R1 或 R2 有多个 primary unique alignment，跳过（不是真正的 unique）
        if len(r1_list) > 1 or len(r2_list) > 1:
            continue
        
        r1 = r1_list[0]
        r2 = r2_list[0]
        
        # 如果指定了 R1 链向，检查是否符合要求
        if r1_strand is not None:
            r1_is_reverse = r1.is_reverse
            if r1_strand == "plus" and r1_is_reverse:
                continue
            if r1_strand == "minus" and not r1_is_reverse:
                continue
        
        # 所有条件都满足，保存这对 alignment
        valid_pairs[read_name] = (r1, r2)
    
    # 写入符合条件的 pair 到输出 BAM
    out_bam = pysam.AlignmentFile(str(output_bam), "wb", template=in_bam)
    
    for read_name, (r1, r2) in valid_pairs.items():
        out_bam.write(r1)
        out_bam.write(r2)
    
    out_bam.close()
    
    # 从 BAM 文件中提取未选中的 reads（包括 unmapped、链向不满足、多重比对的）
    # 重新打开 BAM 文件
    in_bam.close()
    in_bam = pysam.AlignmentFile(str(bam_path), "rb")
    
    # 打开输出 FASTQ 文件（临时文件，排序前）
    output_r1_is_gz = output_r1_fastq.suffix == '.gz'
    output_r2_is_gz = output_r2_fastq.suffix == '.gz'
    
    if output_r1_is_gz:
        out_r1_fq = gzip.open(output_r1_fastq, 'wt')
    else:
        out_r1_fq = open(output_r1_fastq, 'w')
    
    if output_r2_is_gz:
        out_r2_fq = gzip.open(output_r2_fastq, 'wt')
    else:
        out_r2_fq = open(output_r2_fastq, 'w')
    
    # 收集所有未选中的 read names
    valid_read_names = set(valid_pairs.keys())
    
    # 从 BAM 文件中提取未选中的 reads
    # 使用集合记录已处理的 read names，避免重复处理同一个 read 的多个 alignment
    processed_reads = set()
    
    for aln in in_bam:
        read_name = aln.query_name
        
        # 如果这个 read 在 valid_pairs 中，跳过
        if read_name in valid_read_names:
            continue
        
        # 只处理 primary alignment（避免重复处理 secondary 和 supplementary）
        # 但对于 unmapped reads，可能需要处理所有 alignment
        if not aln.is_unmapped and (aln.is_secondary or aln.is_supplementary):
            continue
        
        # 构建唯一标识符（read_name + is_read1/is_read2）
        read_key = (read_name, aln.is_read1, aln.is_read2)
        if read_key in processed_reads:
            continue
        processed_reads.add(read_key)
        
        # 获取序列和质量分数
        seq = aln.query_sequence
        if seq is None or len(seq) == 0:
            continue
        
        # 处理质量分数
        if aln.query_qualities is not None and len(aln.query_qualities) > 0:
            qual = ''.join([chr(q + 33) for q in aln.query_qualities])
        else:
            # 如果没有质量分数，使用 'I' 填充
            qual = 'I' * len(seq)
        
        # 确保质量分数长度与序列长度一致
        if len(qual) != len(seq):
            qual = 'I' * len(seq)
        
        # 处理序列反向互补和质量值反向
        # 如果是 unmapped，直接使用原序列
        # 如果不是 unmapped 且是反向链（is_reverse），需要对序列取反向互补，质量值也要反向
        if not aln.is_unmapped and aln.is_reverse:
            # 序列反向互补
            complement_map = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
            seq = seq.translate(complement_map)[::-1]
            # 质量值反向
            qual = qual[::-1]
        
        # 构建 FASTQ 记录
        header = f"@{read_name}\n"
        fastq_record = f"{header}{seq}\n+\n{qual}\n"
        
        # 根据 flag 判断写入 R1 还是 R2
        if aln.is_read1:
            out_r1_fq.write(fastq_record)
        elif aln.is_read2:
            out_r2_fq.write(fastq_record)
    
    out_r1_fq.close()
    out_r2_fq.close()
    in_bam.close()
    
    # 使用 seqkit pair 重新排序
    logger.info(f"[Seqkit] Pairing and sorting FASTQ files")
    
    # seqkit pair 的输出格式：在输入文件名后面添加 .paired
    # 例如：input.r1.fq -> input.r1.paired.fq
    # 或者：input.r1.fq.gz -> input.r1.paired.fq.gz
    repaired_1 = output_r1_fastq.parent / f"{output_r1_fastq.stem}.paired{output_r1_fastq.suffix}"
    repaired_2 = output_r2_fastq.parent / f"{output_r2_fastq.stem}.paired{output_r2_fastq.suffix}"
    
    # 构建 seqkit pair 命令
    # seqkit pair 会自动在输入文件名后添加 .paired 后缀
    # -o 参数指定输出目录，会在该目录下生成 .paired 文件
    cmd = f"seqkit pair -1 {output_r1_fastq} -2 {output_r2_fastq} -o {output_r1_fastq.parent}/"
    run_cmd(cmd, logger=logger)
    
    # 检查输出文件是否存在（尝试多种可能的命名格式）
    if not repaired_1.exists():
        # 尝试其他可能的命名格式
        # 格式1: {prefix}_1.fq / {prefix}_2.fq
        possible_1 = output_r1_fastq.parent / f"repaired_1{output_r1_fastq.suffix}"
        possible_2 = output_r2_fastq.parent / f"repaired_2{output_r2_fastq.suffix}"
        if possible_1.exists():
            repaired_1 = possible_1
            repaired_2 = possible_2
        else:
            # 格式2: {prefix}.1.fq / {prefix}.2.fq
            possible_1 = output_r1_fastq.parent / f"repaired.1{output_r1_fastq.suffix}"
            possible_2 = output_r2_fastq.parent / f"repaired.2{output_r2_fastq.suffix}"
            if possible_1.exists():
                repaired_1 = possible_1
                repaired_2 = possible_2
            else:
                # 格式3: 基于原始文件名添加 .paired（这是seqkit pair的默认行为）
                # 如果还是没有找到，列出目录中的文件帮助调试
                logger.warning(f"Expected output file not found: {repaired_1}")
                logger.warning(f"Looking for paired files in: {output_r1_fastq.parent}")
                paired_files = list(output_r1_fastq.parent.glob("*.paired.*"))
                if paired_files:
                    logger.warning(f"Found paired files: {paired_files}")
                    # 尝试匹配包含 r1 和 r2 的文件
                    r1_paired = [f for f in paired_files if 'r1' in f.name.lower() or '1' in f.name.lower()]
                    r2_paired = [f for f in paired_files if 'r2' in f.name.lower() or '2' in f.name.lower()]
                    if r1_paired and r2_paired:
                        repaired_1 = r1_paired[0]
                        repaired_2 = r2_paired[0]
    
    if not repaired_1.exists():
        raise FileNotFoundError(
            f"Seqkit pair output file not found: {repaired_1}. "
            f"Expected format: {output_r1_fastq.stem}.paired{output_r1_fastq.suffix}. "
            f"Check seqkit pair output in directory: {output_r1_fastq.parent}"
        )
    if not repaired_2.exists():
        raise FileNotFoundError(
            f"Seqkit pair output file not found: {repaired_2}. "
            f"Expected format: {output_r2_fastq.stem}.paired{output_r2_fastq.suffix}. "
            f"Check seqkit pair output in directory: {output_r2_fastq.parent}"
        )
    
    # 确定最终输出路径（保持原始文件名格式，添加 .repaired 后缀）
    final_repaired_r1 = output_r1_fastq.parent / f"{output_r1_fastq.stem}.repaired{output_r1_fastq.suffix}"
    final_repaired_r2 = output_r2_fastq.parent / f"{output_r2_fastq.stem}.repaired{output_r2_fastq.suffix}"
    
    # 如果输出文件已经是最终路径，不需要移动
    if repaired_1 != final_repaired_r1:
        shutil.move(str(repaired_1), str(final_repaired_r1))
    if repaired_2 != final_repaired_r2:
        shutil.move(str(repaired_2), str(final_repaired_r2))
    
    repaired_r1 = final_repaired_r1
    repaired_r2 = final_repaired_r2
    
    logger.info(f"[Done] Extracted target alignments and dumped others")
    logger.info(f"  - Output BAM: {output_bam}")
    logger.info(f"  - Repaired R1 FASTQ: {repaired_r1}")
    logger.info(f"  - Repaired R2 FASTQ: {repaired_r2}")
    
    return output_bam, repaired_r1, repaired_r2

if __name__ == "__main__":
    bam_path = Path(r"test/test.bam")
    output_bam = Path(r"test/test.target.bam")
    output_r1_fastq = Path(r"test/test.target.r1.fq")
    output_r2_fastq = Path(r"test/test.target.r2.fq")
    r1_strand = "plus"
    extract_target_aln_and_dump_others(bam_path, output_bam, output_r1_fastq, output_r2_fastq, r1_strand)
    print(f"Output BAM: {output_bam}")
    print(f"Output R1 FASTQ: {output_r1_fastq}")
    print(f"Output R2 FASTQ: {output_r2_fastq}")

