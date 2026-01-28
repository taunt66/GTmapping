#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
call_site_FDR.py

从pileup文件中调用修饰位点（m6A/m5C或6mA/5mC），使用二项检验和FDR校正。

支持：
- 碱基类型：A 或 C
- 分子类型：DNA 或 RNA
  - DNA的A是6mA，C是5mC
  - RNA的A是m6A，C是m5C
"""

import logging
from pathlib import Path
from typing import Dict, Tuple, Optional
import scipy.stats
from statsmodels.stats.multitest import multipletests


def collect_background_prob(
    pileup_file: Path,
    base_type: str
) -> Tuple[Dict[str, float], Dict[str, int], Dict[str, int]]:
    """
    统计每条染色体的背景概率（保留率），同时返回统计信息用于输出
    
    参数:
        pileup_file: pileup文件路径
        base_type: 碱基类型，'A' 或 'C'
    
    返回:
        (chr_p, chr_base_total, chr_base_pair_total)
        - chr_p: 每条染色体的背景概率（保留率）
        - chr_base_total: 每条染色体的base计数（A或C）
        - chr_base_pair_total: 每条染色体的base pair计数（AG或CT）
    """
    if base_type.upper() == 'A':
        ref_base = 'A'
        pair_base = 'G'
        chr_base_total = {}
        chr_base_pair_total = {}
    elif base_type.upper() == 'C':
        ref_base = 'C'
        pair_base = 'T'
        chr_base_total = {}
        chr_base_pair_total = {}
    else:
        raise ValueError(f"base_type must be 'A' or 'C', got '{base_type}'")

    with open(pileup_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 5:
                continue
            chr_, pos, ref, strand, base_str = parts
            
            # 跳过表头行：检查pos是否为数字
            try:
                int(pos)  # 检查是否为数字，但不保存
            except ValueError:
                # 如果pos无法转换为整数，可能是表头行，跳过
                continue
            
            bases = base_str.split(",")

            if ref != ref_base:
                continue

            base_count = bases.count(ref_base)
            pair_count = bases.count(pair_base)
            base_pair_total = base_count + pair_count

            if base_pair_total == 0:
                continue

            chr_clean = chr_.split("_AGconvert")[0]
            chr_base_total[chr_clean] = chr_base_total.get(chr_clean, 0) + base_count
            chr_base_pair_total[chr_clean] = chr_base_pair_total.get(chr_clean, 0) + base_pair_total

    # 计算背景概率（保留率）
    chr_p = {}
    for chr_clean in chr_base_total:
        chr_p[chr_clean] = chr_base_total[chr_clean] / chr_base_pair_total[chr_clean]
    
    return chr_p, chr_base_total, chr_base_pair_total


def call_sites_from_pileup(
    pileup_file: Path,
    site_file: Path,
    chr_p: Dict[str, float],
    base_type: str,
    min_pair_cov: int = 10,
    min_base_cov: int = 2,
    min_base_rate: float = 0.1,
    fdr_threshold: float = 0.05,
    logger: Optional[logging.Logger] = None
) -> None:
    """
    从pileup文件中调用修饰位点，使用二项检验和FDR校正
    
    参数:
        pileup_file: pileup文件路径
        site_file: 输出位点文件路径
        chr_p: 每条染色体的背景概率（保留率）
        base_type: 碱基类型，'A' 或 'C'
        min_pair_cov: 最小pair覆盖度（AG或CT）
        min_base_cov: 最小base覆盖度（A或C）
        min_base_rate: 最小base保留率
        fdr_threshold: FDR阈值
        logger: logger对象（可选）
    """
    if logger is None:
        logger = logging.getLogger(__name__)
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('[%(levelname)s] %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
    
    if base_type.upper() == 'A':
        ref_base = 'A'
        pair_base = 'G'
        base_name = 'A'
        pair_name = 'AG'
        rate_name = 'Arate'
    elif base_type.upper() == 'C':
        ref_base = 'C'
        pair_base = 'T'
        base_name = 'C'
        pair_name = 'CT'
        rate_name = 'Crate'
    else:
        raise ValueError(f"base_type must be 'A' or 'C', got '{base_type}'")
    
    logger.info(f"Testing candidate {base_name} sites...")
    results = []  # (chr, pos, strand, base_count, pair_count, rate, pval)

    with open(pileup_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 5:
                continue
            chr_, pos, ref, strand, base_str = parts
            
            # 跳过表头行：检查pos是否为数字
            try:
                pos = int(pos)
            except ValueError:
                # 如果pos无法转换为整数，可能是表头行，跳过
                continue
            
            bases = base_str.split(",")

            if ref != ref_base:
                continue

            base_count = bases.count(ref_base)
            pair_base_count = bases.count(pair_base)
            pair_total = base_count + pair_base_count
            
            if pair_total != 0:
                rate = base_count / pair_total
                if pair_total < min_pair_cov or rate < min_base_rate or base_count < min_base_cov:
                    continue
            else:
                rate = 0

            chr_clean = chr_.split("_AGconvert")[0]
            if chr_clean not in chr_p:
                continue
            p_bg = chr_p[chr_clean]

            # 二项检验：检测该位点保留率是否显著高于背景
            pval = scipy.stats.binom_test(base_count, pair_total, p_bg, alternative="greater")
            results.append((chr_clean, pos, strand, base_count, pair_total, rate, pval))

    # FDR 校正
    if not results:
        logger.warning(f"No candidate {base_name} sites found.")
        return

    pvals = [r[6] for r in results]
    _, padj, _, _ = multipletests(pvals, alpha=fdr_threshold, method="fdr_bh")

    logger.info(f"Writing significant {base_name} sites...")
    with open(site_file, "w") as outf:
        outf.write(f"chr\tpos\tstrand\t{base_name}cov\t{pair_name}cov\t{rate_name}\tPvalue\tPadjust\n")
        for (r, q) in zip(results, padj):
            if q <= fdr_threshold:
                chr_clean, pos, strand, base_count, pair_total, rate, pval = r
                outf.write(
                    f"{chr_clean}\t{pos}\t{strand}\t{base_count}\t{pair_total}\t{rate:.3f}\t{pval:.3e}\t{q:.3e}\n"
                )


def write_chr_conversion_rate(
    chr_p: Dict[str, float],
    output_file: Path,
    base_type: str = "A",
    logger: Optional[logging.Logger] = None
) -> None:
    """
    输出染色体层面的转化率信息（转化率 = 1 - 保留率）
    
    参数:
        chr_p: 每条染色体的背景概率（保留率）
        output_file: 输出文件路径
        base_type: 碱基类型，'A' 或 'C'
        logger: logger对象（可选）
    """
    if logger is None:
        logger = logging.getLogger(__name__)
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('[%(levelname)s] %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
    
    with open(output_file, "w") as outf:
        outf.write("chr\tconversionrate\n")
        for chr_clean in sorted(chr_p.keys()):
            # 转化率 = 1 - 保留率
            conversion_rate = 1 - chr_p[chr_clean]
            outf.write(f"{chr_clean}\t{conversion_rate:.6f}\n")


def get_modification_name(base_type: str, molecule_type: str) -> str:
    """
    根据碱基类型和分子类型获取修饰名称
    
    参数:
        base_type: 碱基类型，'A' 或 'C'
        molecule_type: 分子类型，'DNA' 或 'RNA'
    
    返回:
        修饰名称：'6mA', '5mC', 'm6A', 或 'm5C'
    """
    base_type = base_type.upper()
    molecule_type = molecule_type.upper()
    
    if base_type == 'A':
        if molecule_type == 'DNA':
            return '6mA'
        elif molecule_type == 'RNA':
            return 'm6A'
        else:
            raise ValueError(f"molecule_type must be 'DNA' or 'RNA', got '{molecule_type}'")
    elif base_type == 'C':
        if molecule_type == 'DNA':
            return '5mC'
        elif molecule_type == 'RNA':
            return 'm5C'
        else:
            raise ValueError(f"molecule_type must be 'DNA' or 'RNA', got '{molecule_type}'")
    else:
        raise ValueError(f"base_type must be 'A' or 'C', got '{base_type}'")


def call_site_FDR(
    pileup_file: Path,
    output_dir: Path,
    prefix: str,
    base_type: str,
    molecule_type: str,
    min_pair_cov: int = 10,
    min_base_cov: int = 2,
    min_base_rate: float = 0.1,
    fdr_threshold: float = 0.05,
    logger: Optional[logging.Logger] = None
) -> Tuple[Path, Path]:
    """
    从pileup文件中调用修饰位点，使用二项检验和FDR校正
    
    参数:
        pileup_file: pileup文件路径
        output_dir: 输出目录
        prefix: 输出文件前缀
        base_type: 碱基类型，'A' 或 'C'
        molecule_type: 分子类型，'DNA' 或 'RNA'
          - DNA的A是6mA，C是5mC
          - RNA的A是m6A，C是m5C
        min_pair_cov: 最小pair覆盖度（AG或CT，默认10）
        min_base_cov: 最小base覆盖度（A或C，默认2）
        min_base_rate: 最小base保留率（默认0.1）
        fdr_threshold: FDR阈值（默认0.05）
        logger: logger对象（可选）
    
    返回:
        (site_file_path, conversion_rate_file_path): 输出的位点文件和转化率文件路径的元组
    """
    if logger is None:
        logger = logging.getLogger(__name__)
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('[%(levelname)s] %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
    
    # 验证参数
    base_type = base_type.upper()
    molecule_type = molecule_type.upper()
    
    if base_type not in ['A', 'C']:
        raise ValueError(f"base_type must be 'A' or 'C', got '{base_type}'")
    if molecule_type not in ['DNA', 'RNA']:
        raise ValueError(f"molecule_type must be 'DNA' or 'RNA', got '{molecule_type}'")
    
    # 创建输出目录
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 获取修饰名称
    mod_name = get_modification_name(base_type, molecule_type)
    
    # 定义输出文件路径
    site_file = output_dir / f"{prefix}.{mod_name}.sites.txt"
    conversion_rate_file = output_dir / f"{prefix}.{base_type}.chr_conversion_rate.txt"
    
    # 收集背景概率和统计信息
    logger.info(f"Collecting background probabilities for {base_type}...")
    chr_p, chr_base_total, chr_base_pair_total = collect_background_prob(pileup_file, base_type)
    
    # 输出染色体层面的转化率信息
    logger.info(f"Writing chromosome-level conversion rates for {base_type}...")
    write_chr_conversion_rate(chr_p, conversion_rate_file, base_type, logger)
    
    # 调用修饰位点
    call_sites_from_pileup(
        pileup_file,
        site_file,
        chr_p,
        base_type,
        min_pair_cov,
        min_base_cov,
        min_base_rate,
        fdr_threshold,
        logger
    )
    
    logger.info(f"Analysis completed for {mod_name}!")
    logger.info(f"Output files:")
    logger.info(f"  - {mod_name} sites: {site_file}")
    logger.info(f"  - {base_type} chr conversion rate: {conversion_rate_file}")
    
    return site_file, conversion_rate_file
