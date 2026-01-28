#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_longest_isoform.py

从GTF文件中提取每个基因的最长转录本，并在FASTA文件中提取对应序列。
选择逻辑：对每个基因，选择长度最长的转录本。
"""

from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from typing import Dict, Tuple, List, Optional
import sys
import re


def analyze_gtf_file(gtf_path: Path) -> Tuple[Dict[str, SeqRecord], Dict[str, dict]]:
    """
    解析GTF文件，收集所有基因和转录本信息
    
    参数:
        gtf_path: GTF文件路径
    
    返回:
        gene_infos: {gene_id: SeqRecord}，每个SeqRecord包含转录本features
        transcript_metadata: {transcript_id: {'chr': chr, 'strand': strand, 'gene_id': gene_id, 'gene_name': gene_name}}
    
    注意：
        strand信息用于创建BioPython的FeatureLocation对象（1表示+链，-1表示-链），
        这是正确计算转录本长度所必需的。
    """
    gene_infos = {}
    transcript_infos = {}
    transcript_metadata = {}  # {transcript_id: {'chr': chr, 'strand': strand, 'gene_id': gene_id, 'gene_name': gene_name}}
    gene_metadata = {}  # {gene_id: {'chr': chr, 'strand': strand, 'gene_name': gene_name}}
    
    # 收集所有信息
    with open(gtf_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            chr_name = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes_str = fields[8]
            
            if feature_type == 'gene':
                gene_id_match = re.search(r'gene_id\s+"(.*?)"', attributes_str)
                gene_name_match = re.search(r'gene_name\s+"(.*?)"', attributes_str)
                
                gene_id = gene_id_match.group(1) if gene_id_match else "ELSE"
                gene_name = gene_name_match.group(1) if gene_name_match else gene_id
                
                if gene_id not in gene_infos:
                    gene_infos[gene_id] = SeqRecord(Seq(''), id=gene_id)
                
                # 保存基因元数据
                gene_metadata[gene_id] = {
                    'chr': chr_name,
                    'strand': strand,
                    'gene_name': gene_name
                }
                    
            elif feature_type == 'exon':
                transcript_id_match = re.search(r'transcript_id\s+"(.*?)"', attributes_str)
                gene_id_match = re.search(r'gene_id\s+"(.*?)"', attributes_str)
                gene_name_match = re.search(r'gene_name\s+"(.*?)"', attributes_str)
                
                transcript_id = transcript_id_match.group(1) if transcript_id_match else "ELSE"
                gene_id = gene_id_match.group(1) if gene_id_match else "ELSE"
                gene_name = gene_name_match.group(1) if gene_name_match else gene_id
                
                if 'unknown_transcript' in transcript_id:
                    continue
                
                if transcript_id not in transcript_infos:
                    transcript_infos[transcript_id] = {
                        'location': [(start, end)],
                        'strand': strand,
                        'gene_id': gene_id,
                        'chr': chr_name
                    }
                    # 保存转录本元数据
                    transcript_metadata[transcript_id] = {
                        'chr': chr_name,
                        'strand': strand,
                        'gene_id': gene_id,
                        'gene_name': gene_name
                    }
                else:
                    transcript_infos[transcript_id]['location'].append((start, end))
    
    # 将转录本信息写入基因信息
    for transcript_id, infos in transcript_infos.items():
        exon_locations = [
            FeatureLocation(s, e, strand=1 if infos['strand'] == '+' else -1)
            for s, e in infos['location']
        ]
        
        if len(exon_locations) == 1:
            final_location = exon_locations[0]
        elif len(exon_locations) >= 2:
            final_location = CompoundLocation(exon_locations)
        else:
            continue  # 跳过没有exon的转录本
        
        feature = SeqFeature(location=final_location, id=transcript_id)
        gene_id = infos['gene_id']
        if gene_id in gene_infos:
            gene_infos[gene_id].features.append(feature)
    
    return gene_infos, transcript_metadata


def get_exon_ranges(feature: SeqFeature):
    """
    获取每个exon的起止位置
    
    参数:
        feature: SeqFeature对象
    
    返回:
        生成器，产生(start, end)元组
    """
    if isinstance(feature.location, CompoundLocation):
        for loc in feature.location.parts:
            yield int(loc.start), int(loc.end)
    else:
        yield int(feature.location.start), int(feature.location.end)


def select_longest_transcript(transcripts: Dict[str, Tuple[List[Tuple[int, int]], int]]) -> Tuple[str, List[Tuple[int, int]], int, int, int]:
    """
    选择最长的转录本
    
    参数:
        transcripts: {transcript_id: (exon_list, length)}
    
    返回:
        (transcript_id, exons, start, end, length) 或 None
    """
    if not transcripts:
        return None
    
    # 直接选择最长的转录本
    selected_id = max(transcripts, key=lambda tid: transcripts[tid][1])
    exons, length = transcripts[selected_id]
    start = exons[0][0]
    end = exons[-1][1]
    
    return selected_id, exons, start, end, length


def find_longest_transcript_per_gene(
    gene_infos: Dict[str, SeqRecord],
    transcript_metadata: Dict[str, dict]
) -> Dict[str, dict]:
    """
    找出每个基因的最长转录本
    
    参数:
        gene_infos: {gene_id: SeqRecord}
        transcript_metadata: {transcript_id: {'chr': chr, 'strand': strand, 'gene_id': gene_id, 'gene_name': gene_name}}
    
    返回:
        longest_transcripts: {transcript_id: {
            'gene_id': gene_id,
            'chr': chr,
            'strand': strand,
            'gene_name': gene_name,
            'exons': [(start, end), ...],
            'min_sites': min_start,
            'max_sites': max_end,
            'trans_length': total_length
        }}
    """
    longest_transcripts = {}
    
    for gene_id, record in gene_infos.items():
        transcript_dict = {}  # {transcript_id: (sorted_exons, total_length)}
        
        for feature in record.features:
            transcript_id = feature.id
            
            exons = []
            total_length = 0
            
            for start, end in get_exon_ranges(feature):
                exons.append((start, end))
                total_length += end - start + 1
            
            if total_length > 0:  # 避免空转录本
                # 按起始位置排序
                sorted_exons = sorted(exons, key=lambda x: x[0])
                transcript_dict[transcript_id] = (sorted_exons, total_length)
        
        # 选择最长的转录本
        if transcript_dict:
            result = select_longest_transcript(transcript_dict)
            if result:
                best_transcript_id, exons, min_start, max_end, total_length = result
                
                # 获取转录本元数据
                if best_transcript_id in transcript_metadata:
                    meta = transcript_metadata[best_transcript_id]
                    longest_transcripts[best_transcript_id] = {
                        'gene_id': gene_id,
                        'chr': meta['chr'],
                        'strand': meta['strand'],
                        'gene_name': meta['gene_name'],
                        'exons': exons,
                        'min_sites': min_start,
                        'max_sites': max_end,
                        'trans_length': total_length
                    }
    
    return longest_transcripts


def extract_transcript_id_from_fasta_header(header: str) -> str:
    """
    从FASTA header中提取transcript_id
    
    header格式: >ENST00000832824.1|ENSG00000290825.2|-|-|DDX11L16-260|DDX11L16|1379|lncRNA|
    
    返回:
        transcript_id (如: ENST00000832824.1)
    """
    # 去掉>符号
    header = header.lstrip('>')
    # 按|分割，第一个是transcript_id
    parts = header.split('|')
    if parts:
        return parts[0]
    return ""


def extract_longest_isoform(
    input_fasta: Path, 
    gtf_path: Path, 
    output_fasta: Path,
    output_trans_file: Optional[Path] = None
) -> Tuple[Path, Path]:
    """
    提取每个基因的最长转录本序列，并生成转录本信息文件
    
    参数:
        input_fasta: 输入转录组FASTA文件路径
        gtf_path: GTF注释文件路径
        output_fasta: 输出FASTA文件路径
        output_trans_file: 输出转录本信息文件路径（TSV格式），如果为None则自动生成
    
    返回:
        (output_fasta_path, output_trans_file_path): 输出的FASTA和TSV文件路径的元组
    """
    # 解析GTF文件
    gene_infos, transcript_metadata = analyze_gtf_file(gtf_path)
    
    # 找出每个基因的最长转录本
    longest_transcripts = find_longest_transcript_per_gene(gene_infos, transcript_metadata)
    
    # 创建transcript_id集合（用于快速查找）
    longest_transcript_set = set(longest_transcripts.keys())
    
    # 确保输出目录存在
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    
    # 如果没有指定输出TSV文件路径，自动生成
    if output_trans_file is None:
        output_trans_file = output_fasta.parent / f"{output_fasta.stem}.trans"
    
    output_trans_file.parent.mkdir(parents=True, exist_ok=True)
    
    # 生成转录本信息文件（TSV格式）
    with open(output_trans_file, 'w') as trans_fh:
        # 写入header
        trans_fh.write("Chr\tmin_sites\tmax_sites\texon_starts\texon_ends\tTrans_length\tGene\tStrand\tTranscript_id\n")
        
        # 按transcript_id排序（可选，为了输出的一致性）
        for transcript_id in sorted(longest_transcripts.keys()):
            info = longest_transcripts[transcript_id]
            
            # 格式化外显子信息
            exon_starts = '_'.join([str(e[0]) for e in info['exons']])
            exon_ends = '_'.join([str(e[1]) for e in info['exons']])
            
            # 写入一行
            trans_fh.write(
                f"{info['chr']}\t"
                f"{info['min_sites']}\t"
                f"{info['max_sites']}\t"
                f"{exon_starts}\t"
                f"{exon_ends}\t"
                f"{info['trans_length']}\t"
                f"{info['gene_name']}\t"
                f"{info['strand']}\t"
                f"{transcript_id}\n"
            )
    
    # 读取FASTA文件，提取最长转录本序列
    with open(output_fasta, 'w') as outfh:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            transcript_id = extract_transcript_id_from_fasta_header(record.id)
            
            # 检查这个transcript_id是否是最长转录本
            if transcript_id in longest_transcript_set:
                # 提取序列
                SeqIO.write(record, outfh, 'fasta')
    
    return output_fasta, output_trans_file


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python extract_longest_isoform.py <input_fasta> <gtf_file> <output_fasta> [output_trans_file]", file=sys.stderr)
        sys.exit(1)
    
    input_fasta = Path(sys.argv[1])
    gtf_path = Path(sys.argv[2])
    output_fasta = Path(sys.argv[3])
    output_trans_file = Path(sys.argv[4]) if len(sys.argv) > 4 else None
    
    extract_longest_isoform(input_fasta, gtf_path, output_fasta, output_trans_file)
