#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
生成基于双端测序的配对处理 pileup 结果。

处理逻辑：
- 对于每个位点的每个链向，只统计R1的贡献
- 对于每个R1 alignment，检查其对应的R2是否也在该点相反链向上有覆盖：
  * 如果R2没有覆盖：直接用R1的碱基计数
  * 如果R2也有覆盖：
    - 将R2的碱基互补到当前链向
    - 如果互补后碱基相同：计数
    - 如果不同：取base quality高的
    - 如果base quality相同：不计数（丢弃该pair在该点的贡献）

重要说明：pysam的pileup默认会处理overlapping paired-end reads：
- 默认情况下（ignore_overlaps=True），pysam会检测重叠的read pairs并只保留质量更高的碱基
- 这会导致另一个read的质量值被置为0
- 本脚本已将ignore_overlaps设置为False，以保留原始质量值
- 我们会在后续逻辑中自己处理重叠reads的比较（R1和R2的比较逻辑）

输出格式：chrom  pos  ref_base  strand  bases（1-based 坐标）。
可通过参数指定输出哪些ref碱基（A、G、C、T中的任意1-4个）。

内存优化：每个worker直接写入临时文件，避免内存累积。
"""

import logging
from pathlib import Path
from multiprocessing import Pool
import tempfile
import os
from typing import Optional

import pysam

# 互补碱基映射表（大写），用于快速查找，避免重复调用reverse_complement
_COMP_BASE_MAP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}


def fast_complement(base):
    """快速互补碱基查找，避免字符串操作开销"""
    return _COMP_BASE_MAP.get(base.upper(), 'N')


def pileup_region_paired_to_file(args):
    """
    对单个区域做 pileup，使用配对处理逻辑，直接写入临时文件。
    
    处理逻辑：
    1. 对于每个位点的每个链向，统计主要read（R1或R2，由primary_read参数指定）
    2. 对于每个主要read的alignment，检查其mate pair是否也在该点相反链向上有覆盖：
       - 如果mate没有覆盖：直接用该alignment的碱基计数
       - 如果mate也有覆盖：
         * 将mate的碱基互补到当前链向
         * 如果互补后碱基相同：计数
         * 如果不同：取base quality高的
         * base quality相同：不计数（丢弃该pair在该点的贡献）
    
    参数:
        primary_read: 'R1' 或 'R2'，指定主要统计哪个read
    
    返回临时文件路径
    """
    (
        chrom,
        start,
        end,
        bam_path,
        fasta_path,
        tmp_dir,
        max_depth,
        allowed_ref_bases,
        primary_read,
    ) = args
    
    # 创建临时文件
    tmpfh = None
    tmpname = None
    bamfile = None
    fastafile = None
    
    try:
        tmpfh = tempfile.NamedTemporaryFile(mode='w', delete=False, 
                                            prefix=f"pileup_{chrom}_{start}_", 
                                            suffix=".tmp", 
                                            dir=tmp_dir)
        tmpname = tmpfh.name
        
        bamfile = pysam.AlignmentFile(bam_path, "rb")
        # 使用 pysam.FastaFile 按需读取参考序列，避免内存占用
        fastafile = pysam.FastaFile(str(fasta_path))

        # 注意：pysam 自带的 min_mapq/min_baseq 过滤会丢失低质读，这里先尽量放宽，
        # 在循环内再用用户指定阈值过滤
        try:
            # 关键修复：根据pysam官方文档，ignore_overlaps默认为True，
            # 这会导致pysam检测重叠的read pairs并只保留质量更高的碱基，
            # 将另一个read的质量值置为0。这会导致我们无法获取原始质量值。
            # 因此，我们需要将ignore_overlaps设置为False，以保留原始质量值。
            # 我们会在后续逻辑中自己处理重叠reads的比较。
            pileup_iter = bamfile.pileup(
                chrom,
                start,
                end,
                truncate=True,
                stepper="all",
                max_depth=max_depth,
                flag_filter=0,
                min_base_quality=0,
                min_mapping_quality=0,
                ignore_overlaps=False,  # 设置为False以保留原始质量值
            )
        except TypeError:
            # 对于旧版本的pysam，可能不支持ignore_overlaps参数
            # 尝试使用其他参数
            try:
                pileup_iter = bamfile.pileup(
                    chrom,
                    start,
                    end,
                    truncate=True,
                    stepper="all",
                    max_depth=max_depth,
                    flag_filter=0,
                    min_base_quality=0,
                    min_mapping_quality=0,
                )
            except TypeError:
                # 对于非常旧的版本，使用最简参数
                pileup_iter = bamfile.pileup(
                    chrom,
                    start,
                    end,
                    truncate=True,
                    stepper="all",
                    max_depth=max_depth,
                )

        for pileupcolumn in pileup_iter:
            pos = pileupcolumn.reference_pos  # 0-based
            try:
                # 按需从FASTA文件读取参考序列，避免内存占用
                ref_base = fastafile.fetch(chrom, pos, pos + 1).upper()
                if not ref_base:
                    ref_base = "N"
            except (KeyError, ValueError, IndexError):
                # 如果无法获取参考序列，跳过该位置
                continue
            
            # 第一步：收集该位置所有符合条件的reads信息
            read_infos = []
            _append = read_infos.append  # 缓存方法查找
            
            for pileupread in pileupcolumn.pileups:
                read = pileupread.alignment
                read_name = read.query_name
                is_r1 = read.is_read1
                is_reverse = read.is_reverse
                
                # 检查各种过滤条件
                if pileupread.is_del or pileupread.is_refskip:
                    continue
                
                if read.is_unmapped:
                    continue
                
                if read.is_duplicate:
                    continue
                
                if read.is_qcfail:
                    continue
                
                # 永远不包括secondary比对
                if read.is_secondary:
                    continue
                
                if read.is_supplementary:
                    continue
                
                # 获取query position和碱基
                qp = pileupread.query_position
                if qp is None:
                    continue
                
                try:
                    # 获取原始碱基和质量值
                    # 注意：我们已经将ignore_overlaps设置为False，所以这里获取的是原始质量值
                    query_base = read.query_sequence[qp] if read.query_sequence else None
                    if query_base is None:
                        continue
                    raw_base = query_base.upper()
                    
                    # 获取原始质量值：直接从alignment对象获取
                    # 由于ignore_overlaps=False，这里获取的是原始BAM文件中的质量值
                    if read.query_qualities is not None and qp < len(read.query_qualities):
                        baseq = read.query_qualities[qp]
                    else:
                        baseq = 0
                except (IndexError, TypeError):
                    continue
                
                # 转换到链向方向（使用快速互补查找）
                if read.is_reverse:
                    query_base_on_strand = fast_complement(raw_base)
                    strand = "-"
                else:
                    query_base_on_strand = raw_base
                    strand = "+"
                
                # 收集所有reads（R1和R2）到read_infos，用于后续查找mate进行比较
                _append({
                    'read_name': read_name,
                    'is_r1': is_r1,
                    'base': query_base_on_strand,
                    'baseq': baseq,
                    'strand': strand,
                    'raw_base': raw_base,
                    'is_reverse': is_reverse
                })
            
            # 第二步：建立read_name到reads的映射（用于快速查找mate）
            read_name_to_reads = {}
            strand_to_reads = {'+': [], '-': []}
            
            # 确定主要read类型（用于计数）
            primary_is_r1 = (primary_read.upper() == 'R1')
            
            # 单次循环完成两个索引的构建
            for r in read_infos:
                # 建立read_name+is_r1索引（包含R1和R2，用于查找）
                key = (r['read_name'], r['is_r1'])
                if key not in read_name_to_reads:
                    read_name_to_reads[key] = []
                read_name_to_reads[key].append(r)
                
                # 按链向分组（只添加主要read，用于计数）
                if r['is_r1'] == primary_is_r1:
                    strand_to_reads[r['strand']].append(r)
            
            # 第三步：按链向分组处理
            for strand in ['+', '-']:
                strand_reads = strand_to_reads[strand]
                if not strand_reads:
                    continue
                
                # 第四步：对每个read决定是否计数
                counted_bases = []
                processed_pairs = set()
                _add_base = counted_bases.append  # 缓存方法查找
                
                # 相反链向（用于查找mate）
                opposite_strand = '-' if strand == '+' else '+'
                
                for r in strand_reads:
                    read_name = r['read_name']
                    
                    # strand_to_reads中只包含主要read（R1或R2），所以这里直接处理
                    # 检查是否已经处理过这个pair
                    if read_name in processed_pairs:
                        continue
                    
                    # 查找对应的mate read（read_name相同，is_r1与主要read相反，且在相反链向上）
                    mate_is_r1 = not primary_is_r1
                    mate_key = (read_name, mate_is_r1)
                    mate_reads = read_name_to_reads.get(mate_key, [])
                    # 过滤出相反链向上的mate
                    mate_on_opposite = None
                    for m in mate_reads:
                        if m['strand'] == opposite_strand:
                            mate_on_opposite = m
                            break
                    
                    if mate_on_opposite is None:
                        # mate没有在相反链向上覆盖：直接用主要read的碱基计数
                        _add_base(r['base'])
                        processed_pairs.add(read_name)
                    else:
                        # mate也在相反链向上有覆盖：需要比较
                        # 将mate的碱基互补到当前链向进行比较
                        mate_base_on_current_strand = fast_complement(mate_on_opposite['base'])
                        
                        if r['base'] == mate_base_on_current_strand:
                            # 碱基相同（互补后）：计数
                            _add_base(r['base'])
                        else:
                            # 碱基不同：比较base quality
                            if r['baseq'] > mate_on_opposite['baseq']:
                                _add_base(r['base'])
                            elif mate_on_opposite['baseq'] > r['baseq']:
                                _add_base(mate_base_on_current_strand)
                            # base quality相同：不计数（跳过）
                        
                        processed_pairs.add(read_name)
                
                # 第五步：输出结果（根据allowed_ref_bases过滤），直接写入文件
                if counted_bases:
                    if strand == '+':
                        ref_base_on_strand = ref_base
                    else:  # strand == '-'
                        ref_base_on_strand = fast_complement(ref_base)
                    
                    # 检查ref_base是否在允许的列表中
                    if ref_base_on_strand in allowed_ref_bases:
                        line = f"{chrom}\t{pos + 1}\t{ref_base_on_strand}\t{strand}\t{','.join(counted_bases)}\n"
                        tmpfh.write(line)
                        tmpfh.flush()  # 立即刷新，确保数据写入
            
            # 内存清理：处理完每个位置后立即清理
            del read_infos, read_name_to_reads, strand_to_reads

    except Exception as e:
        # 子进程中无法使用主进程的logger，使用print输出到stderr
        # 这会被主进程的stderr重定向捕获
        import sys
        print(f"ERROR: Error processing chunk {chrom}:{start}-{end}: {e}", file=sys.stderr)
        # 确保临时文件被关闭和删除
        if tmpfh:
            tmpfh.close()
            if tmpname and os.path.exists(tmpname):
                try:
                    os.remove(tmpname)
                except Exception:
                    pass
                tmpname = None
        raise
    finally:
        # 确保资源被清理
        if tmpfh:
            try:
                tmpfh.flush()
                tmpfh.close()
            except Exception:
                pass
        if bamfile:
            try:
                bamfile.close()
            except Exception:
                pass
        if 'fastafile' in locals() and fastafile is not None:
            try:
                fastafile.close()
            except Exception:
                pass
    
    if tmpname is None:
        raise RuntimeError(f"Failed to create temporary file for chunk {chrom}:{start}-{end}")
    
    return tmpname


def get_genome_wide_pileup_paired(
    bam_path: Path,
    fasta_path: Path,
    outputdir: Path,
    thread: int,
    prefix: str,
    tmp_dir: str | None = None,
    chunk_size: int = 500000,
    max_depth: int = 100000,
    allowed_ref_bases: set = None,
    primary_read: str = 'R1',
    logger: Optional[logging.Logger] = None
) -> Path:
    """
    主流程：按区域分块，多进程统计配对处理的 pileup。
    
    处理逻辑：对于每个位点的每个链向，只统计主要read（R1或R2）的贡献。
    当主要read对应的mate在该位置也有覆盖时，使用比较逻辑（比较碱基质量等）。
    
    参数:
        primary_read: 'R1' 或 'R2'，指定主要统计哪个read（默认'R1'）
    
    可指定输出哪些ref碱基（A、G、C、T中的任意1-4个）。
    
    内存优化：每个worker直接写入临时文件，最后合并所有临时文件。
    """
    if allowed_ref_bases is None:
        allowed_ref_bases = {'A'}  # 默认只输出A
    
    # 验证primary_read参数
    primary_read_upper = primary_read.upper()
    if primary_read_upper not in ['R1', 'R2']:
        raise ValueError(f"primary_read must be 'R1' or 'R2', got '{primary_read}'")
    primary_read = primary_read_upper
    
    # 如果没有传入logger，使用默认的logging模块
    if logger is None:
        logger = logging.getLogger("pipeline")
        if not logger.handlers:
            # 如果没有handler，创建一个简单的handler
            handler = logging.StreamHandler()
            formatter = logging.Formatter('[%(levelname)s] %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
    
    outputdir.mkdir(parents=True, exist_ok=True)
    
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    chrom_lengths = dict(zip(bamfile.references, bamfile.lengths))
    bamfile.close()

    # 处理全部基因组
    regions = [(chrom, 0, length) for chrom, length in chrom_lengths.items()]

    # 准备临时目录
    if tmp_dir is None:
        tmp_dir = tempfile.gettempdir()
    os.makedirs(tmp_dir, exist_ok=True)

    # 构建任务列表
    tasks = []
    for chrom, length_start, length_end in regions:
        cur = length_start
        while cur < length_end:
            end = min(cur + chunk_size, length_end)
            tasks.append(
                (
                    chrom,
                    cur,
                    end,
                    str(bam_path),
                    str(fasta_path),
                    tmp_dir,
                    max_depth,
                    allowed_ref_bases,
                    primary_read,
                )
            )
            cur = end

    logger.info(
        f"Processing paired-end pileup using {thread} threads "
        f"(chunk_size={chunk_size}, max_depth={max_depth}, primary_read={primary_read})"
    )
    logger.info(f"Total tasks: {len(tasks)}")
    
    # 多进程处理，使用 imap 以更好地控制内存和任务分配
    tmp_paths = []
    with Pool(processes=thread) as pool:
        # 使用 imap 而不是 map，可以更好地控制内存使用
        # chunksize=1 确保任务逐个分配，避免一次性加载所有任务到内存
        results = pool.imap(pileup_region_paired_to_file, tasks, chunksize=1)
        # 处理返回结果
        for i, result in enumerate(results):
            tmp_paths.append(result)
            # 每处理100个任务输出一次进度
            if (i + 1) % 100 == 0:
                logger.info(f"Completed {i + 1}/{len(tasks)} tasks")
    
    # 合并所有临时文件
    output_file = outputdir / f"{prefix}.pileup"
    logger.info(f"Merging {len(tmp_paths)} temporary files into {output_file}")
    
    try:
        with open(output_file, "w") as outfh:
            # 写入表头
            outfh.write("chrom\tpos\tref_base\tstrand\tbases\n")
            
            # 按顺序合并临时文件
            for i, tmp_path in enumerate(tmp_paths):
                if not os.path.exists(tmp_path):
                    logger.warning(f"Temporary file {tmp_path} does not exist, skipping")
                    continue
                
                try:
                    with open(tmp_path, "r") as tmpfh:
                        # 使用块读取以减少内存占用
                        while True:
                            chunk = tmpfh.read(8192)  # 8KB chunks
                            if not chunk:
                                break
                            outfh.write(chunk)
                except Exception as e:
                    logger.error(f"Error reading tmp file {tmp_path}: {e}")
                    continue
                finally:
                    # 删除临时文件
                    try:
                        os.remove(tmp_path)
                    except Exception as e:
                        logger.warning(f"Could not remove tmp file {tmp_path}: {e}")
                
                # 进度日志
                if (i + 1) % 100 == 0:
                    logger.info(f"Merged {i + 1}/{len(tmp_paths)} temporary files")
    
    except Exception as e:
        logger.error(f"Error writing output file {output_file}: {e}")
        # 清理剩余的临时文件
        for tmp_path in tmp_paths:
            try:
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)
            except Exception:
                pass
        raise

    logger.info("Paired-end pileup writing completed.")
    return output_file
