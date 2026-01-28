# 快速开始指南

本文档提供快速上手指南，帮助用户快速运行流程。

## 前置要求

1. **安装必需工具**：
   - STAR (>=2.7.0)
   - Bowtie2 (>=2.4.0)
   - samtools (>=1.10)
   - seqkit (>=2.0.0)

2. **安装Python依赖**：
   ```bash
   pip install -r requirements.txt
   ```

3. **准备参考文件**：
   - 基因组FASTA文件
   - 转录组FASTA文件（RNA流程需要）
   - GTF注释文件（RNA流程需要）

## 快速开始

### 1. 下载参考文件（可选）

使用内置下载工具下载参考文件：

```bash
# 下载E.coli参考文件
python -m final_PE_pipeline.scripts.download \
    --tooldir /path/to/tool/dir \
    --species ecoli

# 下载人类参考文件
python -m final_PE_pipeline.scripts.download \
    --tooldir /path/to/tool/dir \
    --species human
```

### 2. 运行流程

#### DNA流程示例

```bash
python -m final_PE_pipeline.run \
    --r1 R1.fastq.gz \
    --r2 R2.fastq.gz \
    --genome-fa /path/to/genome.fa \
    --trans-fa /path/to/transcriptome.fa \
    --gtf /path/to/annotation.gtf \
    --tooldir /path/to/tool/dir \
    --outputdir /path/to/output \
    --prefix sample_dna \
    --threads 24 \
    --type DNA
```

#### RNA流程示例

```bash
python -m final_PE_pipeline.run \
    --r1 R1.fastq.gz \
    --r2 R2.fastq.gz \
    --genome-fa /path/to/genome.fa \
    --trans-fa /path/to/transcriptome.fa \
    --gtf /path/to/annotation.gtf \
    --tooldir /path/to/tool/dir \
    --outputdir /path/to/output \
    --prefix sample_rna \
    --threads 24 \
    --type RNA
```

### 3. 使用Dry Run模式测试

在正式运行前，建议先使用dry run模式检查流程：

```bash
python -m final_PE_pipeline.run \
    ... [其他参数] ... \
    --dryrun
```

## 参数说明

### 必需参数

- `--r1`: R1 FASTQ文件路径
- `--r2`: R2 FASTQ文件路径
- `--genome-fa`: 基因组FASTA文件路径
- `--trans-fa`: 转录组FASTA文件路径（RNA流程必需）
- `--gtf`: GTF注释文件路径（RNA流程必需）
- `--tooldir`: 工具目录（存储索引和中间文件）
- `--outputdir`: 输出目录
- `--prefix`: 输出文件前缀
- `--type`: 分子类型（`DNA` 或 `RNA`）

### 可选参数

- `--threads`: 线程数（默认：8）
- `--dryrun`: 干运行模式
- `--filter-target-bases`: 过滤目标碱基（默认：DNA为`T G`，RNA为`A C`）
- `--filter-cutoffs`: 过滤阈值（默认：`3`）
- `--min-pair-cov`: 最小pair覆盖度（默认：10）
- `--min-base-cov`: 最小碱基覆盖度（默认：2）
- `--min-base-rate`: 最小碱基比率（默认：0.1）
- `--fdr-threshold`: FDR阈值（默认：0.05）

## 输出文件

运行完成后，在`outputdir`中会生成：

- `{prefix}.pileup`: Pileup文件
- `{prefix}.{modification}.sites.txt`: 修饰位点文件
  - DNA: `{prefix}.6mA.sites.txt`, `{prefix}.5mC.sites.txt`
  - RNA: `{prefix}.m6A.sites.txt`, `{prefix}.m5C.sites.txt`
- `{prefix}.{base}.chr_conversion_rate.txt`: 染色体转化率文件

## 常见问题

### Q: 如何知道流程是否正常运行？

A: 检查日志文件（位于`tooldir/logs/`），所有操作都会记录在日志中。

### Q: DNA流程需要转录组文件吗？

A: 不需要。DNA流程只需要基因组FASTA文件。但命令行参数仍需要提供`--trans-fa`和`--gtf`（可以为空文件），这些参数在DNA流程中不会被使用。

### Q: 如何调整过滤参数？

A: 使用`--filter-target-bases`和`--filter-cutoffs`参数。例如：
```bash
--filter-target-bases A C \
--filter-cutoffs 3 3
```

### Q: 流程运行时间大概多久？

A: 取决于数据量和硬件配置。一般：
- 小数据集（<10M reads）：几小时
- 中等数据集（10-100M reads）：半天到一天
- 大数据集（>100M reads）：1-3天

### Q: 如何并行运行多个样本？

A: 每个样本独立运行，可以并行执行。注意：
- 确保`tooldir`中的索引文件可以被多个进程共享
- 每个样本使用不同的`--outputdir`和`--prefix`

## 下一步

- 查看 [README.md](README.md) 了解详细功能
- 查看 [INTERFACE.md](INTERFACE.md) 了解接口规范
- 查看 [CODE_STYLE.md](CODE_STYLE.md) 了解代码风格
- 查看 [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md) 了解项目结构

