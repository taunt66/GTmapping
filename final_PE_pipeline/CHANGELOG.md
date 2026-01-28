# 更新日志

## [2.1.0] - 接口统一和文档更新

### 重大变更
- **文件位置调整**：
  - `scripts/pre_process.py` → `pre_process.py`（移至根目录）
  - `scripts/download.py` → `download.py`（移至根目录）
- **接口统一**：所有生成BAM的函数现在都自动进行排序和索引
- **文档全面更新**：更新所有相关文档以反映最新接口和文件结构

### 新增
- **自动BAM处理**：
  - `run_cmd_and_logging.py`: 新增`sort_and_index_bam()`函数
  - 所有mapping函数自动对输出的BAM进行排序和索引
  - `transTogenome.py`: 转换后的BAM自动排序和索引
  - `merge_bam.py`: 合并后的BAM自动排序和索引
  - `return_origin_bam_and_filter.py`: 过滤后的BAM自动排序和索引
  - `extract_target_aln_dump_others.py`: 提取后的BAM自动建立索引

- **pre_process.py增强**：
  - RNA流程：新增完整的预处理流程（与new_library_preprocess.py一致）
    - Step 1: 筛选R2的11-13bp是AAA的pair，提取UMI，去AAA接头
    - Step 2: seqkit rmdup -s去重
    - Step 3: cutadapt -u 23 -m 20切除UMI和额外trim
  - 固定参数：RNA流程的UMI长度（10bp）、trim长度（3bp）、接头序列（AAA）不再暴露为命令行参数

- **return_origin_bam_and_filter.py增强**：
  - 修复质量值丢失问题：在恢复序列时保留原始质量值
  - 支持多碱基类型和阈值的统一过滤（`target_bases`和`cutoffs`为列表）

### 修改
- **所有BAM生成函数**：
  - `mapping.py`: STAR和bowtie2输出的BAM自动排序和索引
  - `transTogenome.py`: 转换后的BAM自动排序和索引
  - `merge_bam.py`: 合并后的BAM自动排序和索引
  - `return_origin_bam_and_filter.py`: 过滤后的BAM自动排序和索引
  - `extract_target_aln_dump_others.py`: 提取后的BAM自动建立索引

- **文档更新**：
  - `INTERFACE.md`: 完整更新所有模块的接口列表，包括新增的`sort_and_index_bam()`函数
  - `PROJECT_STRUCTURE.md`: 反映文件位置变化，更新所有模块说明
  - `README.md`: 更新目录结构和工作流程说明
  - `CHANGELOG.md`: 记录本次统一接口的变更

- **接口统一**：
  - 所有函数都遵循统一的参数顺序：必需参数在前，可选参数在后，`logger`参数在最后
  - 所有函数都包含完整的类型提示和文档字符串
  - 所有命令执行都通过`run_cmd()`函数

### 修复
- **return_origin_bam_and_filter.py**: 修复质量值丢失问题（pysam在修改序列时会清空质量值）
- **extract_target_aln_dump_others.py**: 改进`seqkit pair`输出文件检测逻辑
- **Spike_in_analysis.py**: 修复BAM索引问题，确保pileup前BAM已排序和索引

### 技术改进
- 统一的BAM处理流程：所有BAM文件在生成后都会自动排序和索引
- 改进的错误处理：所有函数都包含输入验证和输出验证
- 改进的日志记录：所有操作都记录到统一的日志系统

## [2.0.0] - DNA/RNA流程分离和代码风格统一

### 重大变更
- **DNA/RNA流程分离**：根据`--type`参数自动选择不同的处理流程
  - DNA流程：两步基因组mapping，过滤T/G，按R2计数
  - RNA流程：三步mapping（两次基因组+一次转录组），过滤A/C，按R1计数
- **代码风格统一**：创建`CODE_STYLE.md`，统一所有脚本的代码风格

### 新增
- `CODE_STYLE.md`：详细的代码风格规范文档
- DNA流程支持：独立的文件准备和mapping流程
- `pre_process.py`：FASTQ预处理模块（DNA/RNA特定处理）
- `download.py`：参考文件下载工具

### 修改
- **run.py**:
  - 分离`check_and_prepare_files_dna()`和`check_and_prepare_files_rna()`
  - 根据type选择不同的mapping流程
  - 根据type设置不同的过滤参数和pileup参数
- **README.md**: 更新为反映DNA/RNA流程差异
- **INTERFACE.md**: 添加新模块接口说明

### 参数变更
- `--filter-target-bases`: 默认值根据type自动设置（DNA: T G，RNA: A C）
- `--filter-cutoffs`: 默认值为3（对所有碱基）

## [1.0.0] - 整理为GitHub版本

### 新增
- 添加了`README.md`文档：包含项目介绍、安装说明、使用方法等
- 添加了`requirements.txt`：Python依赖包列表
- 添加了`.gitignore`：Git忽略文件配置
- 添加了`INTERFACE.md`：接口规范文档
- 添加了`CHANGELOG.md`：更新日志
- 添加了`__init__.py`：模块初始化文件

### 修改
- **merge_bam.py**: 
  - 修改为使用`samtools`命令（通过`run_cmd()`）替代`pysam.merge/sort/index`
  - 添加了`logger`参数支持
  - 所有命令执行都通过`run_cmd_and_logging`模块进行日志记录

- **run.py**: 
  - 更新`merge_bams`调用，添加`logger`参数
  - 所有流程步骤都已通过统一接口调用

### 统一接口
- 所有脚本函数现在都支持`logger`参数
- 所有Linux命令都通过`run_cmd()`执行，确保日志记录统一
- 统一的参数命名和类型提示
- 详细的文档字符串（中文）

### 代码规范
- 统一使用类型提示
- 统一的错误处理
- 统一的日志记录方式

## 技术栈

- **Python**: 3.8+
- **主要工具**: STAR, Bowtie2, samtools, seqkit, cutadapt
- **Python库**: pysam, biopython, pandas, scipy, statsmodels

## 注意事项

- 所有命令执行必须通过`run_cmd()`函数，以确保日志记录
- 所有函数都应接受`logger`参数（可选），如果没有提供则使用默认logger
- 所有生成BAM的函数都会自动进行排序和索引，无需手动处理
- 大型数据文件（FASTA、GTF、BAM等）不应提交到Git仓库（已在.gitignore中配置）
