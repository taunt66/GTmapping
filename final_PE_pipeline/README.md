# GT-Mapping Pipeline

双端测序数据处理流程，用于检测DNA/RNA修饰位点（6mA/5mC或m6A/m5C）。

**快速开始**: 查看 [QUICK_START.md](QUICK_START.md) 获取快速上手指南。

## 功能概述

本流程实现了从原始双端FASTQ文件到修饰位点检测的完整处理流程，支持**DNA**和**RNA**两种分子类型：

### DNA流程
1. **FASTQ预处理**：R1去除5'端poly-G（>=7个G），去除接头
2. **文件准备**：碱基转换（R1: T→C, G→A; R2: A→G, C→T）、构建索引（仅基因组）
3. **两次Mapping**：STAR基因组TC_GA（正链）、STAR基因组AG_CT（负链）
4. **BAM合并**：合并两个mapping结果
5. **原始BAM还原**：使用索引文件还原原始序列并过滤（T和G）
6. **Pileup分析**：按R2计数，生成A和C的覆盖度统计
7. **位点检测**：检测6mA和5mC修饰位点

### RNA流程
1. **FASTQ预处理**：筛选R2的11-13bp是AAA的pair，提取UMI，去接头，去重，trim
2. **文件准备**：提取最长转录本、碱基转换（R1: A→G, C→T; R2: T→C, G→A）、构建索引
3. **三次Mapping**：STAR基因组AG_CT（正链）、STAR基因组TC_GA（负链）、bowtie2转录组AG_CT（正链）
4. **坐标转换**：转录组坐标转基因组坐标
5. **BAM合并**：合并三个mapping结果
6. **原始BAM还原**：使用索引文件还原原始序列并过滤（A和C）
7. **Pileup分析**：按R1计数，生成A和C的覆盖度统计
8. **位点检测**：检测m6A和m5C修饰位点

## 目录结构

详细的项目结构说明请参考 [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md)。

```
final_PE_pipeline/
├── README.md                    # 本文件（使用说明）
├── CODE_STYLE.md                # 代码风格规范
├── INTERFACE.md                 # 接口规范文档
├── CHANGELOG.md                 # 更新日志
├── PROJECT_STRUCTURE.md         # 项目结构说明
├── QUICK_START.md               # 快速开始指南
├── WORKFLOW_LOGIC.md            # 工作流程逻辑详解
├── requirements.txt             # Python依赖
├── .gitignore                   # Git忽略文件
├── download_url.txt             # 下载链接（参考数据）
├── run.py                       # 主流程入口脚本
├── pre_process.py               # FASTQ预处理（DNA/RNA）
├── download.py                  # 参考文件下载工具
├── Spike_in_analysis.py         # Spike-in分析脚本
├── test_spikein.fa              # Spike-in测试FASTA文件
├── test_spikein.index           # Spike-in测试索引文件
└── scripts/                     # 功能脚本目录
    ├── run_cmd_and_logging.py   # 命令执行和日志记录（核心）
    ├── build_index.py           # 构建索引（STAR/bowtie2）
    ├── mapping.py               # Mapping（STAR/bowtie2）
    ├── convert_fastq_fasta.py  # 碱基转换（FASTQ/FASTA）
    ├── extract_longest_isoform.py  # 提取最长转录本（RNA）
    ├── extract_target_aln_dump_others.py  # 提取目标比对
    ├── transTogenome.py         # 转录组坐标转基因组坐标（RNA）
    ├── merge_bam.py             # BAM文件合并
    ├── return_origin_bam_and_filter.py  # 还原原始BAM并过滤
    ├── pileup_PE.py             # 双端Pileup分析
    ├── call_site_FDR.py         # 修饰位点检测（FDR校正）
    └── EXTRACT_FILES_EXPLANATION.md  # 提取文件说明文档
```

## 安装要求

### 系统要求

- Linux/macOS
- Python 3.8+

### 必需工具

- **STAR** (>=2.7.0): 用于RNA/DNA测序数据比对
- **Bowtie2** (>=2.4.0): 用于转录组比对
- **samtools** (>=1.10): 用于BAM文件处理
- **seqkit** (>=2.0.0): 用于FASTQ文件处理

### Python依赖

```bash
pip install -r requirements.txt
```

## 使用方法

**快速开始**: 查看 [QUICK_START.md](QUICK_START.md) 获取详细的使用示例和常见问题解答。

### 基本用法

```bash
python -m final_PE_pipeline.run \
    --r1 R1.fastq.gz \
    --r2 R2.fastq.gz \
    --genome-fa genome.fa \
    --trans-fa transcriptome.fa \
    --gtf annotation.gtf \
    --tooldir /path/to/tool/dir \
    --outputdir /path/to/output \
    --prefix sample \
    --threads 24 \
    --type RNA
```

### 参数说明

- `--r1`: R1 FASTQ文件路径（必需）
- `--r2`: R2 FASTQ文件路径（必需）
- `--genome-fa`: 基因组FASTA文件路径（必需）
- `--trans-fa`: 转录组FASTA文件路径（必需）
- `--gtf`: GTF注释文件路径（必需）
- `--tooldir`: 工具目录（项目根目录，用于存储索引和中间文件）（必需）
- `--outputdir`: 输出目录（必需）
- `--prefix`: 输出文件前缀（必需）
- `--threads`: 线程数（默认：8）
- `--type`: 分子类型，`DNA` 或 `RNA`（必需）
- `--dryrun`: 干运行模式，仅打印操作不执行（可选）
- `--filter-target-bases`: 过滤目标碱基列表（默认：DNA为`T G`，RNA为`A C`）
- `--filter-cutoffs`: 过滤阈值列表（默认：`3`，对所有碱基）
- `--min-pair-cov`: 最小pair覆盖度（默认：10）
- `--min-base-cov`: 最小碱基覆盖度（默认：2）
- `--min-base-rate`: 最小碱基比率（默认：0.1）
- `--fdr-threshold`: FDR阈值（默认：0.05）

### 干运行模式（Dry Run）

用于预览流程而不实际执行：

```bash
python -m final_PE_pipeline.run \
    ... [其他参数] ... \
    --dryrun
```

## 工作流程详解

### DNA流程

#### Step 0: FASTQ预处理
- **DNA流程**：R1去除5'端poly-G（>=7个G），去除接头
- **RNA流程**：筛选R2的11-13bp是AAA的pair，提取UMI，去接头，去重，trim

#### Step 1: FASTQ转换和索引生成
- **R1转换**：T→C, G→A
- **R2转换**：A→G, C→T
- 生成索引文件用于后续还原

#### Step 2: 文件准备
- **碱基转换**：
  - 基因组TC_GA转换（T→C, G→A）
  - 基因组AG_CT转换（A→G, C→T）
- **构建索引**：
  - STAR索引（基因组TC_GA、AG_CT）

#### Step 3: Mapping流程
1. **STAR基因组TC_GA（正链）**：提取R1映射到正链的pair
2. **STAR基因组AG_CT（负链）**：剩余FASTQ映射，提取R1映射到负链的pair
3. **BAM合并**：合并两个mapping结果

#### Step 3: 原始BAM还原和过滤
- 使用索引文件将BAM中的碱基还原为原始序列
- 过滤含有过量T或G的pair（默认cutoff=3）

#### Step 4: Pileup分析
- 按R2计数
- 统计A和C的覆盖度信息

#### Step 5: 位点检测
- 检测6mA和5mC修饰位点
- 使用二项检验和FDR校正

### RNA流程

#### Step 0: FASTQ预处理
- 筛选R2的11-13bp是AAA的pair
- 提取R1和R2各10bp UMI拼成20bp
- 对R1剩余部分去除AAA接头
- seqkit rmdup -s去重
- cutadapt -u 23 -m 20切除UMI和额外trim

#### Step 1: FASTQ转换和索引生成
- **R1转换**：A→G, C→T
- **R2转换**：T→C, G→A
- 生成索引文件用于后续还原

#### Step 2: 文件准备
- **提取最长转录本**：从GTF文件中提取每个基因的最长转录本
- **碱基转换**：
  - 基因组AG_CT转换（A→G, C→T）
  - 基因组TC_GA转换（T→C, G→A）
  - 转录组AG_CT转换（A→G, C→T）
- **构建索引**：
  - STAR索引（基因组AG_CT、TC_GA）
  - Bowtie2索引（转录组AG_CT）

#### Step 3: Mapping流程
1. **STAR基因组AG_CT（正链）**：提取R1映射到正链的pair
2. **STAR基因组TC_GA（负链）**：剩余FASTQ映射，提取R1映射到负链的pair
3. **Bowtie2转录组AG_CT（正链）**：剩余FASTQ映射，提取R1映射到正链的pair
4. **坐标转换**：将转录组坐标的BAM转换为基因组坐标
5. **BAM合并**：合并三个mapping结果

#### Step 4: 原始BAM还原和过滤
- 使用索引文件将BAM中的碱基还原为原始序列
- 过滤含有过量A或C的pair（默认cutoff=3）

#### Step 5: Pileup分析
- 按R1计数
- 统计A和C的覆盖度信息

#### Step 6: 位点检测
- 检测m6A和m5C修饰位点
- 使用二项检验和FDR校正（Benjamini-Hochberg方法）

## 输出文件

### 主要输出

- `{prefix}.pileup`: Pileup文件（位点覆盖度统计）
- `{prefix}.{modification}.sites.txt`: 修饰位点文件（如`{prefix}.m6A.sites.txt`）
- `{prefix}.{base}.chr_conversion_rate.txt`: 染色体转化率文件

### 中间文件

- `alignment/`: BAM文件和中间比对结果
- `fastq/`: 转换后的FASTQ和索引文件
- `reference/`: 转换后的参考序列
- `annotations/`: 转录本信息文件
- `index/`: 比对索引文件

## 脚本模块说明

### run_cmd_and_logging.py

统一的命令执行和日志记录模块：
- `setup_project_logger()`: 设置项目级logger
- `run_cmd()`: 执行shell命令并记录日志
- `get_logger()`: 获取logger实例

**所有Linux命令必须通过`run_cmd()`执行，以确保日志记录。**

### 其他脚本

所有脚本都提供函数接口，支持：
- `logger`参数（可选）：用于日志记录
- 统一的参数命名和类型
- 详细的文档字符串
- 自动BAM排序和索引（所有生成BAM的函数）

### 主要脚本模块

- **pre_process.py**（根目录）：FASTQ预处理，支持DNA和RNA两种流程
- **download.py**（根目录）：参考文件下载工具，支持E.coli和人类参考文件
- **Spike_in_analysis.py**（根目录）：Spike-in序列甲基化水平分析

## 开发说明

### 代码规范

详细的代码风格规范请参考 [CODE_STYLE.md](CODE_STYLE.md)。

主要规范：
1. **统一接口**：所有脚本函数接受`logger`参数
2. **命令执行**：所有Linux命令通过`run_cmd()`执行
3. **类型提示**：使用Python类型提示
4. **文档字符串**：所有函数都有详细的中文文档
5. **导入顺序**：标准库 → 第三方库 → 本地模块

### 测试

建议使用`--dryrun`参数测试流程：

```bash
python -m final_PE_pipeline.run --dryrun ...
```

## 许可证

[根据项目需要添加许可证信息]

## 贡献

欢迎提交Issue和Pull Request。

## 联系方式

[根据项目需要添加联系方式]

