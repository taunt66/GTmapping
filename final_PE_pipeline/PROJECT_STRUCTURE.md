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
├── QUICK_START.md               # 快速开始指南
├── WORKFLOW_LOGIC.md            # 工作流程逻辑详解
├── requirements.txt             # Python依赖包列表
├── .gitignore                   # Git忽略文件配置
├── download_url.txt             # 参考数据下载链接
├── __init__.py                  # Python包初始化文件
├── run.py                       # 主流程入口脚本
├── pre_process.py               # FASTQ预处理脚本（DNA/RNA）
├── download.py                  # 参考文件下载工具
├── Spike_in_analysis.py        # Spike-in分析脚本
├── test_spikein.fa             # Spike-in测试FASTA文件
├── test_spikein.index          # Spike-in测试索引文件
└── scripts/                     # 功能脚本目录
    ├── __init__.py              # 脚本模块初始化
    ├── run_cmd_and_logging.py   # 命令执行和日志记录（核心模块）
    ├── build_index.py           # 索引构建（STAR/bowtie2）
    ├── mapping.py               # 比对执行（STAR/bowtie2）
    ├── convert_fastq_fasta.py   # 碱基转换（FASTQ/FASTA）
    ├── extract_longest_isoform.py  # 提取最长转录本（RNA流程）
    ├── extract_target_aln_dump_others.py  # 提取目标比对并输出未比对reads
    ├── transTogenome.py         # 转录组坐标转基因组坐标（RNA流程）
    ├── merge_bam.py             # BAM文件合并
    ├── return_origin_bam_and_filter.py  # 还原原始BAM并过滤
    ├── pileup_PE.py             # 双端Pileup分析
    ├── call_site_FDR.py         # 修饰位点检测（FDR校正）
    └── EXTRACT_FILES_EXPLANATION.md  # 提取文件说明文档
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

- **QUICK_START.md**: 快速开始指南，包含使用示例和常见问题

- **WORKFLOW_LOGIC.md**: 工作流程逻辑详解，包含DNA和RNA流程的详细步骤

### 主流程脚本

- **run.py**: 主流程入口脚本
  - 解析命令行参数
  - 根据type（DNA/RNA）选择不同流程
  - 协调各模块执行
  - 支持dryrun模式

- **pre_process.py**: FASTQ预处理脚本（位于根目录）
  - DNA流程：R1去除5'端poly-G，去除接头
  - RNA流程：UMI提取、去接头、去重、trim
  - 自动识别gzip压缩格式
  - 支持命令行接口和函数接口

- **download.py**: 参考文件下载工具（位于根目录）
  - 下载E.coli和人类参考文件
  - 初始化项目目录结构
  - 支持断点续传和进度显示

- **Spike_in_analysis.py**: Spike-in分析脚本（位于根目录）
  - 分析spike-in序列的甲基化水平
  - 自动转换FASTA和FASTQ
  - 执行mapping和pileup分析

### 核心模块（scripts/）

#### run_cmd_and_logging.py
**功能**: 统一的命令执行和日志记录
- `setup_project_logger()`: 设置项目级logger，捕获stdout/stderr
- `run_cmd()`: 执行shell命令并记录日志
- `run_cmd_list()`: 顺序执行多个命令
- `sort_and_index_bam()`: 对BAM文件进行排序并建立索引
- `get_logger()`: 获取logger实例
- **重要性**: 所有Linux命令必须通过此模块执行

#### build_index.py
**功能**: 构建比对索引
- `build_index()`: 统一接口（STAR/bowtie2）
- `build_star_index()`: 构建STAR索引（支持GTF文件）
- `build_bowtie2_index()`: 构建bowtie2索引
- `calculate_genome_length()`: 计算基因组长度（用于STAR参数）

#### mapping.py
**功能**: 执行比对
- `mapping()`: 统一接口（STAR/bowtie2）
- `mapping_star()`: STAR比对（支持RNA/DNA模式，自动检测gzip）
- `mapping_bowtie2()`: Bowtie2比对
- **自动功能**: 所有mapping输出的BAM都会自动排序和索引

#### convert_fastq_fasta.py
**功能**: 碱基转换
- `change_fastq()`: 转换FASTQ并生成索引文件
- `change_fasta()`: 转换FASTA文件
- 支持多种转换规则（A→G, C→T, T→C, G→A等）
- 自动识别gzip压缩格式

#### extract_longest_isoform.py
**功能**: 提取最长转录本（RNA流程）
- `extract_longest_isoform()`: 从GTF提取最长转录本
- `analyze_gtf_file()`: 解析GTF文件
- `find_longest_transcript_per_gene()`: 找出每个基因的最长转录本
- 生成转录本信息文件（TSV格式）

#### extract_target_aln_dump_others.py
**功能**: 提取目标比对并输出未比对reads
- `extract_target_aln_and_dump_others()`: 提取符合条件的pair，输出未选中的reads
- `is_primary_unique()`: 判断是否为primary unique比对
- 支持按R1链向过滤（plus/minus）
- 对反向链reads进行序列反向互补和质量值反向
- 使用seqkit pair重新配对FASTQ
- **自动功能**: 输出的BAM会自动建立索引

#### transTogenome.py
**功能**: 转录组坐标转基因组坐标（RNA流程）
- `change_trans_to_genome()`: 转换BAM坐标
- `translate_cigar()`: 转换CIGAR字符串
- `change_segment()`: 转换单个read segment
- 处理双端数据，正确设置mate信息
- **自动功能**: 输出的BAM会自动排序和索引

#### merge_bam.py
**功能**: BAM文件合并
- `merge_bams()`: 合并多个BAM文件
- `read_headers()`: 读取BAM header
- `process_and_write_bam()`: 处理并写入BAM
- 统一染色体名称（去掉后缀如`_AGconvert`）
- 使用samtools命令执行
- **自动功能**: 合并后的BAM会自动排序和索引

#### return_origin_bam_and_filter.py
**功能**: 还原原始BAM并过滤
- `filter_and_return_origin_bam()`: 还原原始序列并过滤
- `parse_index_file()`: 解析索引文件
- `restore_sequence()`: 还原序列（保留原始质量值）
- `count_base_in_read()`: 统计碱基数量
- 支持多碱基类型和阈值的统一过滤
- **自动功能**: 过滤后的BAM会自动排序和索引

#### pileup_PE.py
**功能**: 双端Pileup分析
- `get_genome_wide_pileup_paired()`: 全基因组pileup
- `pileup_region_paired_to_file()`: 区域pileup（多进程）
- 支持指定主要read（R1/R2）
- 支持指定允许的参考碱基（A/C/G/T）
- 处理paired-end重叠reads
- 使用多进程处理，内存优化

#### call_site_FDR.py
**功能**: 修饰位点检测（FDR校正）
- `call_site_FDR()`: 主函数，调用修饰位点
- `get_modification_name()`: 获取修饰名称（6mA/5mC/m6A/m5C）
- `collect_background_prob()`: 收集背景概率
- `call_sites_from_pileup()`: 从pileup调用位点
- 使用二项检验和Benjamini-Hochberg FDR校正

## 工作流程依赖关系

### DNA流程
```
run.py
  ├── pre_process.py (FASTQ预处理)
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
  ├── pre_process.py (FASTQ预处理)
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
- `pre_process.py` → `run_cmd_and_logging.py` (命令执行)
- `download.py` → `run_cmd_and_logging.py` (日志记录)

## 数据流

### DNA流程数据流
```
原始FASTQ (R1, R2)
  ↓ [pre_process]
预处理后FASTQ
  ↓ [convert_fastq_fasta]
转换后FASTQ + 索引文件
  ↓ [mapping × 2]
两个BAM文件（已排序和索引）
  ↓ [merge_bam]
合并BAM（已排序和索引）
  ↓ [return_origin_bam_and_filter]
还原并过滤的BAM（已排序和索引）
  ↓ [pileup_PE]
Pileup文件
  ↓ [call_site_FDR]
修饰位点文件
```

### RNA流程数据流
```
原始FASTQ (R1, R2)
  ↓ [pre_process]
预处理后FASTQ
  ↓ [convert_fastq_fasta]
转换后FASTQ + 索引文件
  ↓ [extract_longest_isoform]
最长转录本FASTA + 转录本信息文件
  ↓ [mapping × 3]
三个BAM文件（已排序和索引）
  ↓ [transTogenome]
转录组BAM转基因组坐标（已排序和索引）
  ↓ [merge_bam]
合并BAM（已排序和索引）
  ↓ [return_origin_bam_and_filter]
还原并过滤的BAM（已排序和索引）
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
│   ├── *.AG_CTconvert.fa
│   ├── *.TC_GAconvert.fa
│   └── *.longest.fa (RNA流程)
├── annotations/        # 转录本信息文件（RNA流程）
│   └── *.trans
├── index/              # 比对索引文件
│   ├── *_STAR_index/
│   └── *_bowtie2_index.*
└── logs/               # 日志文件
    ├── pipeline.log.{timestamp}
    └── preprocess_*.log.{timestamp}

outputdir/
├── alignment/          # BAM文件和中间比对结果
│   ├── *.bam
│   ├── *.bam.bai
│   └── *.remaining_*.paired.fq.gz
├── fastq/              # 转换后的FASTQ和索引文件
│   ├── *_converted.fq.gz
│   └── *_index_*.txt
├── {prefix}.pileup     # Pileup文件
├── {prefix}.{mod}.sites.txt  # 修饰位点文件
└── {prefix}.{base}.chr_conversion_rate.txt  # 转化率文件
```

## 自动功能

以下操作会自动执行，无需手动调用：

1. **BAM排序和索引**：
   - `mapping.py`: 所有mapping输出的BAM自动排序和索引
   - `transTogenome.py`: 转换后的BAM自动排序和索引
   - `merge_bam.py`: 合并后的BAM自动排序和索引
   - `return_origin_bam_and_filter.py`: 过滤后的BAM自动排序和索引
   - `extract_target_aln_dump_others.py`: 提取后的BAM自动建立索引

2. **文件格式识别**：
   - `pre_process.py`: 自动识别gzip压缩格式
   - `convert_fastq_fasta.py`: 自动识别gzip压缩格式
   - `mapping.py`: 自动检测gzip并添加`--readFilesCommand zcat`

3. **临时文件清理**：
   - 所有模块都会自动清理临时文件
   - 中间文件会在处理完成后删除
