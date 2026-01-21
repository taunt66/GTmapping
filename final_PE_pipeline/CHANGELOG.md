# 更新日志

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
- **主要工具**: STAR, Bowtie2, samtools, seqkit
- **Python库**: pysam, biopython, pandas, scipy, statsmodels

## 注意事项

- 所有命令执行必须通过`run_cmd()`函数，以确保日志记录
- 所有函数都应接受`logger`参数（可选），如果没有提供则使用默认logger
- 大型数据文件（FASTA、GTF、BAM等）不应提交到Git仓库（已在.gitignore中配置）

