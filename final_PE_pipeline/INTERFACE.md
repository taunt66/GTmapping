# 接口规范文档

本文档定义了`final_PE_pipeline`中所有脚本模块的统一接口规范。

**注意**：详细的代码风格规范请参考 [CODE_STYLE.md](CODE_STYLE.md)。

## 基本原则

1. **所有Linux命令必须通过`run_cmd()`执行**：确保命令执行和日志记录的统一
2. **所有函数接受`logger`参数**：用于日志记录，类型为`Optional[logging.Logger]`，默认值为`None`
3. **使用类型提示**：所有函数参数和返回值都应包含类型提示
4. **文档字符串**：所有公开函数都应包含详细的中文文档字符串

## 函数签名规范

### 基本格式

```python
def function_name(
    required_param: Type,
    optional_param: Type = default_value,
    logger: Optional[logging.Logger] = None
) -> ReturnType:
    """
    函数功能描述
    
    参数:
        required_param: 必需参数描述
        optional_param: 可选参数描述（默认值：xxx）
        logger: logger对象（可选），如果不提供则使用默认logger
    
    返回:
        返回值描述
    """
    if logger is None:
        logger = get_logger()  # 从run_cmd_and_logging导入
    # 函数实现...
```

### Logger处理

所有函数都应支持`logger`参数，如果没有提供，则使用默认logger：

```python
from .run_cmd_and_logging import get_logger

def my_function(param: str, logger: Optional[logging.Logger] = None):
    if logger is None:
        logger = get_logger()
    logger.info("...")
```

### 命令执行

所有Linux命令必须通过`run_cmd()`执行：

```python
from .run_cmd_and_logging import run_cmd

# 正确：通过run_cmd执行
run_cmd("samtools sort -o output.bam input.bam", logger=logger)

# 错误：直接使用subprocess或os.system
# subprocess.run("samtools sort ...")  # 不允许
# os.system("samtools sort ...")  # 不允许
# pysam.merge(...)  # 对于命令操作，应该使用samtools命令
```

## 模块接口列表

### run_cmd_and_logging.py

**核心命令执行模块**

- `setup_project_logger(log_dir: Path, log_name: str, level: int) -> logging.Logger`
  - 设置项目级logger
- `run_cmd(cmd: str, logger: Optional[logging.Logger] = None, ...) -> subprocess.CompletedProcess`
  - 执行shell命令并记录日志
- `get_logger() -> logging.Logger`
  - 获取默认logger实例

### build_index.py

**索引构建模块**

- `build_index(tool: str, genome_fasta: Path, output_path: Path, threads: int = 8, gtf_file: Optional[Path] = None, sjdb_overhang: int = 100, logger: Optional[logging.Logger] = None) -> Path`
  - 统一接口：构建STAR或bowtie2索引
- `build_star_index(genome_fasta: Path, gtf_file: Path, output_dir: Path, threads: int = 8, sjdb_overhang: int = 100, logger: Optional[logging.Logger] = None) -> Path`
  - 构建STAR索引
- `build_bowtie2_index(genome_fasta: Path, output_prefix: Path, threads: int = 1, logger: Optional[logging.Logger] = None) -> Path`
  - 构建bowtie2索引

### mapping.py

**比对模块**

- `mapping(tool: str, r1_fastq: Path, r2_fastq: Path, index_path: Path, output_path: Path, type: str = "RNA", threads: int = 8, output_prefix: Optional[Path] = None, mismatch_ratio: float = 0.03, logger: Optional[logging.Logger] = None) -> Path`
  - 统一接口：执行STAR或bowtie2比对
- `mapping_star(...) -> Path`
  - STAR比对
- `mapping_bowtie2(...) -> Path`
  - Bowtie2比对

### extract_target_aln_dump_others.py

**比对提取和转储模块**

- `extract_target_aln_and_dump_others(bam_path: Path, output_bam: Path, output_r1_fastq: Path, output_r2_fastq: Path, r1_strand: Optional[str] = None, logger: Optional[logging.Logger] = None) -> Tuple[Path, Path, Path]`
  - 提取目标比对并输出未比对reads

### merge_bam.py

**BAM合并模块**

- `merge_bams(bams: List[Path], outputdir: Path, prefix: str, logger: Optional[logging.Logger] = None) -> Path`
  - 合并多个BAM文件（使用samtools命令）

### transTogenome.py

**坐标转换模块**

- `change_trans_to_genome(trans_bam: Path, target_genome_fa: Path, longest_transfile: Path, outputbam: Path) -> Path`
  - 将转录组坐标BAM转换为基因组坐标BAM

### extract_longest_isoform.py

**转录本提取模块**

- `extract_longest_isoform(input_fasta: Path, gtf_path: Path, output_fasta: Path, output_trans_file: Optional[Path] = None) -> Tuple[Path, Path]`
  - 提取最长转录本

### convert_fastq_fasta.py

**碱基转换模块**

- `change_fastq(input_fastq: Path, output_dir: Path, convert_list: List[str]) -> Tuple[Path, Path]`
  - 转换FASTQ文件并生成索引
- `change_fasta(input_fasta: Path, output_fasta: Path, convert_list: List[str]) -> Path`
  - 转换FASTA文件

### return_origin_bam_and_filter.py

**BAM还原和过滤模块**

- `filter_and_return_origin_bam(bam_path: Path, r1_index_file: Path, r2_index_file: Path, origin_bam: Path, filter_bam: Path, target_base: str = 'A', cutoff: int = 3, logger: Optional[logging.Logger] = None) -> Tuple[Path, Path]`
  - 还原原始BAM并过滤

### pileup_PE.py

**Pileup分析模块**

- `get_genome_wide_pileup_paired(bam_path: Path, fasta_path: Path, outputdir: Path, thread: int, prefix: str, tmp_dir: str | None = None, chunk_size: int = 500000, max_depth: int = 100000, allowed_ref_bases: set = None, primary_read: str = 'R1') -> Path`
  - 生成双端配对处理的pileup
  - `primary_read`: 'R1' 或 'R2'，指定主要统计哪个read（DNA流程使用R2，RNA流程使用R1）

### pre_process.py

**FASTQ预处理模块**

- `pre_process_fastq(r1_in: Path, r2_in: Path, r1_out: Path, r2_out: Path, type: str, threads: int = 4, logger: Optional[logging.Logger] = None)`
  - 根据分子类型（DNA/RNA）对FASTQ进行预处理
  - DNA: R1去除5'端9个G，R2去除3'端9个C
  - RNA: 直接复制文件

### download.py

**参考文件下载模块**

- `download_ecoli_references(tooldir: Path, logger: Optional[logging.Logger] = None)`
  - 下载E.coli K12 MG1655参考文件
- `download_human_references(tooldir: Path, logger: Optional[logging.Logger] = None)`
  - 下载人类GRCh38参考文件
- `initialize_project_directories(tooldir: Path, logger: Optional[logging.Logger] = None)`
  - 初始化项目目录结构

### call_site_FDR.py

**位点检测模块**

- `call_site_FDR(pileup_file: Path, output_dir: Path, prefix: str, base_type: str, molecule_type: str, min_pair_cov: int = 10, min_base_cov: int = 2, min_base_rate: float = 0.1, fdr_threshold: float = 0.05, logger: Optional[logging.Logger] = None) -> Tuple[Path, Path]`
  - 从pileup调用修饰位点（FDR校正）

## 命名规范

1. **函数名**：使用小写字母和下划线（snake_case）
2. **参数名**：使用小写字母和下划线（snake_case）
3. **类型参数**：使用大写字母开头的驼峰命名（PascalCase）

## 导入规范

模块间导入应优先使用相对导入：

```python
from .run_cmd_and_logging import run_cmd, get_logger
```

如果相对导入失败，可以回退到绝对导入或fallback：

```python
try:
    from .run_cmd_and_logging import run_cmd, get_logger
except ImportError:
    from run_cmd_and_logging import run_cmd, get_logger
```

## 错误处理

所有函数应该：
1. 验证输入参数
2. 提供清晰的错误消息
3. 记录关键操作到日志

## 示例

### 标准函数实现

```python
from pathlib import Path
from typing import Optional
import logging
from .run_cmd_and_logging import run_cmd, get_logger

def process_file(
    input_file: Path,
    output_file: Path,
    threads: int = 8,
    logger: Optional[logging.Logger] = None
) -> Path:
    """
    处理文件的示例函数
    
    参数:
        input_file: 输入文件路径
        output_file: 输出文件路径
        threads: 线程数（默认：8）
        logger: logger对象（可选）
    
    返回:
        输出文件路径
    """
    if logger is None:
        logger = get_logger()
    
    # 验证输入
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    # 确保输出目录存在
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # 执行命令（通过run_cmd）
    cmd = f"tool --input {input_file} --output {output_file} --threads {threads}"
    run_cmd(cmd, logger=logger)
    
    # 验证输出
    if not output_file.exists():
        raise FileNotFoundError(f"Output file was not created: {output_file}")
    
    logger.info(f"Processing completed: {output_file}")
    return output_file
```

