# 项目结构说明

本文档详细说明`final_PE_pipeline`项目的目录结构和各文件功能。

## 目录结构

```
final_PE_pipeline/
├── README.md                    # 项目主文档（使用说明、安装、流程介绍）
├── CODE_STYLE.md                # 代码风格规范
├── INTERFACE.md                 # 接口规范文档
├── CHANGELOG.md                 # 更新日志
├── PROJECT_STRUCTURE.md         # 本文件（项目结构说明）
├── requirements.txt             # Python依赖包列表
├── .gitignore                   # Git忽略文件配置
├── __init__.py                  # Python包初始化文件
├── run.py                       # 主流程入口脚本
├── download_url.txt             # 参考数据下载链接
├── Spike_in_analysis.py         # Spike-in分析脚本（待实现）
└── scripts/                     # 功能脚本目录
    ├── __init__.py              # 脚本模块初始化
    ├── run_cmd_and_logging.py   # 命令执行和日志记录（核心模块）
    ├── build_index.py           # 索引构建（STAR/bowtie2）
    ├── mapping.py               # 比对执行（STAR/bowtie2）
    ├── convert_fastq_fasta.py   # 碱基转换（FASTQ/FASTA）
    ├── pre_process.py           # FASTQ预处理（DNA/RNA特定处理）
    ├── extract_longest_isoform.py  # 提取最长转录本（RNA流程）
    ├── extract_target_aln_dump_others.py  # 提取目标比对并输出未比对reads
    ├── transTogenome.py         # 转录组坐标转基因组坐标（RNA流程）
    ├── merge_bam.py             # BAM文件合并
    ├── return_origin_bam_and_filter.py  # 还原原始BAM并过滤
    ├── pileup_PE.py             # 双端Pileup分析
    ├── call_site_FDR.py         # 修饰位点检测（FDR校正）
    └── download.py              # 参考文件下载工具
```

## 文件功能说明

### 核心文档

- **README.md**: 项目主文档，包含：
  - 功能概述
  - 安装要求
  - 使用方法
  - 工作流程详解（DNA/RNA流程）
  - 输出文件说明

- **CODE_STYLE.md**: 代码风格规范，包含：
  - 文件头格式
  - 导入顺序规范
  - 函数文档字符串格式
  - 命名规范
  - 错误处理规范

- **INTERFACE.md**: 接口规范文档，包含：
  - 函数签名规范
  - Logger处理规范
  - 命令执行规范
  - 各模块接口列表

- **CHANGELOG.md**: 更新日志，记录版本变更

- **PROJECT_STRUCTURE.md**: 本文件，项目结构说明

### 主流程脚本

- **run.py**: 主流程入口脚本
  - 解析命令行参数
  - 根据type（DNA/RNA）选择不同流程
  - 协调各模块执行
  - 支持dryrun模式

### 核心模块（scripts/）

#### run_cmd_and_logging.py
**功能**: 统一的命令执行和日志记录
- `setup_project_logger()`: 设置项目级logger
- `run_cmd()`: 执行shell命令并记录日志
- `get_logger()`: 获取logger实例
- **重要性**: 所有Linux命令必须通过此模块执行

#### build_index.py
**功能**: 构建比对索引
- `build_index()`: 统一接口（STAR/bowtie2）
- `build_star_index()`: 构建STAR索引
- `build_bowtie2_index()`: 构建bowtie2索引
- `calculate_genome_length()`: 计算基因组长度（用于STAR参数）

#### mapping.py
**功能**: 执行比对
- `mapping()`: 统一接口（STAR/bowtie2）
- `mapping_star()`: STAR比对（支持RNA/DNA模式）
- `mapping_bowtie2()`: Bowtie2比对

#### convert_fastq_fasta.py
**功能**: 碱基转换
- `change_fastq()`: 转换FASTQ并生成索引文件
- `change_fasta()`: 转换FASTA文件
- 支持多种转换规则（A→G, C→T等）

#### pre_process.py
**功能**: FASTQ预处理（DNA/RNA特定处理）
- `pre_process_fastq()`: 根据type进行预处理
  - DNA: R1去除5'端9个G，R2去除3'端9个C
  - RNA: 直接复制文件

#### extract_longest_isoform.py
**功能**: 提取最长转录本（RNA流程）
- `extract_longest_isoform()`: 从GTF提取最长转录本
- `analyze_gtf_file()`: 解析GTF文件
- `find_longest_transcript_per_gene()`: 找出每个基因的最长转录本

#### extract_target_aln_dump_others.py
**功能**: 提取目标比对并输出未比对reads
- `extract_target_aln_and_dump_others()`: 提取符合条件的pair，输出未选中的reads
- `is_primary_unique()`: 判断是否为primary unique比对
- 支持按R1链向过滤

#### transTogenome.py
**功能**: 转录组坐标转基因组坐标（RNA流程）
- `change_trans_to_genome()`: 转换BAM坐标
- `translate_cigar()`: 转换CIGAR字符串
- `change_segment()`: 转换单个read segment

#### merge_bam.py
**功能**: BAM文件合并
- `merge_bams()`: 合并多个BAM文件
- `read_headers()`: 读取BAM header
- `process_and_write_bam()`: 处理并写入BAM
- 使用samtools命令执行

#### return_origin_bam_and_filter.py
**功能**: 还原原始BAM并过滤
- `filter_and_return_origin_bam()`: 还原原始序列并过滤
- `parse_index_file()`: 解析索引文件
- `restore_sequence()`: 还原序列
- `count_base_in_read()`: 统计碱基数量

#### pileup_PE.py
**功能**: 双端Pileup分析
- `get_genome_wide_pileup_paired()`: 全基因组pileup
- `pileup_region_paired_to_file()`: 区域pileup（多进程）
- 支持指定主要read（R1/R2）
- 处理paired-end重叠reads

#### call_site_FDR.py
**功能**: 修饰位点检测（FDR校正）
- `call_site_FDR()`: 主函数，调用修饰位点
- `get_modification_name()`: 获取修饰名称（6mA/5mC/m6A/m5C）
- `_collect_background_prob()`: 收集背景概率
- `_call_sites_from_pileup()`: 从pileup调用位点
- 使用二项检验和FDR校正

#### download.py
**功能**: 参考文件下载工具
- `download_ecoli_references()`: 下载E.coli参考文件
- `download_human_references()`: 下载人类参考文件
- `initialize_project_directories()`: 初始化项目目录结构
- `download_file()`: 通用下载函数（带进度显示）

## 工作流程依赖关系

### DNA流程
```
run.py
  ├── convert_fastq_fasta.py (FASTQ转换)
  ├── build_index.py (构建索引)
  ├── mapping.py (两次mapping)
  ├── merge_bam.py (合并BAM)
  ├── return_origin_bam_and_filter.py (还原和过滤)
  ├── pileup_PE.py (Pileup分析)
  └── call_site_FDR.py (位点检测)
```

### RNA流程
```
run.py
  ├── convert_fastq_fasta.py (FASTQ转换)
  ├── extract_longest_isoform.py (提取最长转录本)
  ├── build_index.py (构建索引)
  ├── mapping.py (三次mapping)
  ├── transTogenome.py (坐标转换)
  ├── merge_bam.py (合并BAM)
  ├── return_origin_bam_and_filter.py (还原和过滤)
  ├── pileup_PE.py (Pileup分析)
  └── call_site_FDR.py (位点检测)
```

## 模块依赖关系

所有模块都依赖：
- `run_cmd_and_logging.py`: 命令执行和日志记录

其他依赖：
- `mapping.py` → `build_index.py` (需要索引)
- `transTogenome.py` → `extract_longest_isoform.py` (需要转录本信息)
- `return_origin_bam_and_filter.py` → `convert_fastq_fasta.py` (需要索引文件)

## 数据流

### DNA流程数据流
```
原始FASTQ (R1, R2)
  ↓ [convert_fastq_fasta]
转换后FASTQ + 索引文件
  ↓ [mapping × 2]
两个BAM文件
  ↓ [merge_bam]
合并BAM
  ↓ [return_origin_bam_and_filter]
还原并过滤的BAM
  ↓ [pileup_PE]
Pileup文件
  ↓ [call_site_FDR]
修饰位点文件
```

### RNA流程数据流
```
原始FASTQ (R1, R2)
  ↓ [convert_fastq_fasta]
转换后FASTQ + 索引文件
  ↓ [mapping × 3]
三个BAM文件（两个基因组 + 一个转录组）
  ↓ [transTogenome]
转录组BAM转基因组坐标
  ↓ [merge_bam]
合并BAM
  ↓ [return_origin_bam_and_filter]
还原并过滤的BAM
  ↓ [pileup_PE]
Pileup文件
  ↓ [call_site_FDR]
修饰位点文件
```

## 输出目录结构

运行后，`tooldir`和`outputdir`会包含以下结构：

```
tooldir/
├── reference/          # 转换后的参考序列
├── annotations/        # 转录本信息文件（RNA流程）
├── index/              # 比对索引文件
└── logs/               # 日志文件

outputdir/
├── alignment/          # BAM文件和中间比对结果
├── fastq/              # 转换后的FASTQ和索引文件
├── {prefix}.pileup     # Pileup文件
├── {prefix}.{mod}.sites.txt  # 修饰位点文件
└── {prefix}.{base}.chr_conversion_rate.txt  # 转化率文件
```

