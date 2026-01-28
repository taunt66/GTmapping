# 接口规范文档

本文档定义了`final_PE_pipeline`中所有脚本模块的统一接口规范。

**注意**：详细的代码风格规范请参考 [CODE_STYLE.md](CODE_STYLE.md)。

## 基本原则

1. **所有Linux命令必须通过`run_cmd()`执行**：确保命令执行和日志记录的统一
2. **所有函数接受`logger`参数**：用于日志记录，类型为`Optional[logging.Logger]`，默认值为`None`
3. **使用类型提示**：所有函数参数和返回值都应包含类型提示
4. **文档字符串**：所有公开函数都应包含详细的中文文档字符串
5. **参数顺序**：必需参数在前，可选参数在后，`logger`参数放在最后

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
```

## 模块接口列表

### run_cmd_and_logging.py

**核心命令执行模块**（位于`scripts/`目录）

- `setup_project_logger(log_dir: Path, log_name: str = "pipeline.log", level: int = logging.INFO) -> logging.Logger`
  - 设置项目级logger，捕获stdout和stderr
  - 自动添加时间戳到日志文件名

- `run_cmd(cmd: str, logger: Optional[logging.Logger] = None, shell: bool = True, check: bool = True, cwd: Optional[Path] = None, env: Optional[dict] = None) -> subprocess.CompletedProcess`
  - 执行shell命令并记录日志
  - 自动记录命令、输出和返回码

- `run_cmd_list(cmd_list: List[str], logger: Optional[logging.Logger] = None, ...) -> List[subprocess.CompletedProcess]`
  - 顺序执行多个命令

- `sort_and_index_bam(bam_path: Path, logger: Optional[logging.Logger] = None, already_sorted: bool = False, threads: Optional[int] = None) -> Path`
  - 对BAM文件进行排序并建立索引
  - 如果`already_sorted=True`，只建立索引

- `get_logger() -> logging.Logger`
  - 获取默认logger实例

### pre_process.py

**FASTQ预处理模块**（位于根目录）

- `preprocess_paired_fastq(r1_in: Path, r2_in: Path, r1_out: Path, r2_out: Path, type: str, threads: int = 4, logger: Optional[logging.Logger] = None, tooldir: Optional[Path] = None) -> None`
  - 根据分子类型（DNA/RNA）对FASTQ进行预处理
  - **DNA流程**：
    - Step 1: 切掉R1 5'端的所有G（>=7个G才保留），R2保持不变
    - Step 2: 使用cutadapt去除接头（R1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC, R2: CCCCCCC）
  - **RNA流程**：
    - Step 1: 筛选R2的11-13bp是AAA的pair，提取R1和R2各10bp UMI拼成20bp，对R1剩余部分去除AAA接头
    - Step 2: 使用seqkit rmdup -s去重
    - Step 3: 使用cutadapt -u 23 -m 20切除UMI和额外trim
  - 自动识别gzip压缩格式

- `trim_polyg_and_filter_pairs(r1_in: Path, r2_in: Path, r1_out: Path, r2_out: Path, min_g_count: int = 7, logger: Optional[logging.Logger] = None) -> Tuple[int, int]`
  - DNA专用：切掉R1 5'端的所有G并过滤pair

- `remove_adapters_with_cutadapt(r1_in: Path, r2_in: Path, r1_out: Path, r2_out: Path, threads: int = 24, logger: Optional[logging.Logger] = None) -> None`
  - 使用cutadapt去除接头序列

### download.py

**参考文件下载模块**（位于根目录）

- `download_file(url: str, output_path: Path, logger: Optional[logging.Logger] = None) -> Path`
  - 下载文件（支持断点续传，显示进度）

- `decompress_file(compressed_path: Path, output_path: Optional[Path] = None, logger: Optional[logging.Logger] = None) -> Path`
  - 解压.gz文件

- `initialize_project_directories(tooldir: Path, logger: Optional[logging.Logger] = None) -> None`
  - 初始化项目目录结构（reference/, annotations/, index/, logs/等）

- `download_ecoli_references(tooldir: Path, logger: Optional[logging.Logger] = None) -> Tuple[Path, Path]`
  - 下载E.coli K12 MG1655参考文件（基因组FASTA和GTF）
  - 返回: (genome_fa, gtf_file)

- `download_human_references(tooldir: Path, logger: Optional[logging.Logger] = None) -> Tuple[Path, Path, Path]`
  - 下载人类GRCh38参考文件（基因组FASTA、转录组FASTA和GTF）
  - 返回: (genome_fa, trans_fa, gtf_file)

### build_index.py

**索引构建模块**（位于`scripts/`目录）

- `build_index(tool: str, genome_fasta: Path, output_path: Path, threads: int = 8, gtf_file: Optional[Path] = None, sjdb_overhang: int = 100, logger: Optional[logging.Logger] = None) -> Path`
  - 统一接口：构建STAR或bowtie2索引
  - `tool`: 'STAR' 或 'bowtie2'
  - 对于STAR，如果提供`gtf_file`，会使用GTF文件构建splicing-aware索引

- `build_star_index(genome_fasta: Path, gtf_file: Optional[Path], output_dir: Path, threads: int = 8, sjdb_overhang: int = 100, logger: Optional[logging.Logger] = None) -> Path`
  - 构建STAR索引
  - 自动计算`genomeSAindexNbases`参数

- `build_bowtie2_index(genome_fasta: Path, output_prefix: Path, threads: int = 1, logger: Optional[logging.Logger] = None) -> Path`
  - 构建bowtie2索引

### mapping.py

**比对模块**（位于`scripts/`目录）

- `mapping(tool: str, r1_fastq: Path, r2_fastq: Path, index_path: Path, output_path: Path, type: str = "RNA", threads: int = 8, output_prefix: Optional[Path] = None, mismatch_ratio: float = 0.03, logger: Optional[logging.Logger] = None) -> Path`
  - 统一接口：执行STAR或bowtie2比对
  - 自动对输出的BAM进行排序和索引
  - 对于STAR，`output_path`作为`output_prefix`使用
  - 对于bowtie2，`output_path`作为输出SAM路径（会自动转换为BAM）

- `mapping_star(r1_fastq: Path, r2_fastq: Path, index_dir: Path, output_prefix: Path, type: str = "RNA", threads: int = 24, mismatch_ratio: float = 0.03, logger: Optional[logging.Logger] = None) -> Path`
  - STAR比对
  - 自动检测gzip压缩格式并添加`--readFilesCommand zcat`
  - RNA模式：使用splicing-aware参数
  - DNA模式：使用no-splice参数（禁用splicing）
  - 输出已排序的BAM，并自动建立索引

- `mapping_bowtie2(r1_fastq: Path, r2_fastq: Path, index_prefix: Path, output_sam: Path, type: str = "RNA", threads: int = 8, logger: Optional[logging.Logger] = None) -> Path`
  - Bowtie2比对
  - 自动将SAM转换为BAM，删除SAM文件
  - 对输出的BAM进行排序并建立索引

### convert_fastq_fasta.py

**碱基转换模块**（位于`scripts/`目录）

- `change_fastq(input_fastq: Path, output_dir: Path, convert_list: List[str], output_fastq: Optional[Path] = None) -> Tuple[Path, Path]`
  - 转换FASTQ文件并生成索引文件
  - `convert_list`: 转换规则列表，如`['A:G', 'C:T']`
  - 返回: (转换后FASTQ路径, 索引文件路径)
  - 自动识别gzip压缩格式
  - 索引文件格式：`read_id\tA_index\tC_index`（用`_`分隔多个位置）

- `change_fasta(input_fasta: Path, output_fasta: Path, convert_list: List[str]) -> Path`
  - 转换FASTA文件
  - 返回: 转换后FASTA文件路径

### extract_longest_isoform.py

**转录本提取模块**（位于`scripts/`目录）

- `extract_longest_isoform(input_fasta: Path, gtf_path: Path, output_fasta: Path, output_trans_file: Optional[Path] = None) -> Tuple[Path, Path]`
  - 提取每个基因的最长转录本序列
  - 生成转录本信息文件（TSV格式）
  - 返回: (输出FASTA路径, 转录本信息文件路径)

### extract_target_aln_dump_others.py

**比对提取和转储模块**（位于`scripts/`目录）

- `extract_target_aln_and_dump_others(bam_path: Path, output_bam: Path, output_r1_fastq: Path, output_r2_fastq: Path, r1_strand: Optional[str] = None, logger: Optional[logging.Logger] = None) -> Tuple[Path, Optional[Path], Optional[Path]]`
  - 提取目标比对并输出未比对reads
  - `r1_strand`: 'plus'、'minus'或None（提取所有）
  - 只提取unique primary alignment（NH:i:1）
  - 对非unmapped且反向的reads进行序列反向互补和质量值反向
  - 使用seqkit pair重新配对FASTQ
  - 自动对输出的BAM建立索引
  - 如果输出FASTQ路径为`/dev/null`，则跳过FASTQ提取
  - 返回: (输出BAM路径, R1 FASTQ路径, R2 FASTQ路径)，如果跳过FASTQ提取，后两个为None

### transTogenome.py

**坐标转换模块**（位于`scripts/`目录）

- `change_trans_to_genome(trans_bam: Path, target_genome_fa: Path, longest_transfile: Path, outputbam: Path) -> Path`
  - 将转录组坐标BAM转换为基因组坐标BAM
  - 处理双端数据，正确设置mate信息
  - 自动对输出的BAM进行排序并建立索引
  - 返回: 输出BAM文件路径

### merge_bam.py

**BAM合并模块**（位于`scripts/`目录）

- `merge_bams(bams: List[Path], outputdir: Path, prefix: str, logger: Optional[logging.Logger] = None) -> Path`
  - 合并多个BAM文件
  - 统一染色体名称（去掉后缀如`_AGconvert`）
  - 使用samtools命令执行
  - 自动对合并后的BAM进行排序并建立索引
  - 返回: 合并并排序后的BAM文件路径

### return_origin_bam_and_filter.py

**BAM还原和过滤模块**（位于`scripts/`目录）

- `filter_and_return_origin_bam(bam_path: Path, r1_index_file: Path, r2_index_file: Path, origin_bam: Path, filter_bam: Path, target_bases: List[str] = None, cutoffs: List[int] = None, logger: Optional[logging.Logger] = None) -> Tuple[Path, Path]`
  - 还原原始BAM并过滤
  - `target_bases`: 要检查的碱基列表（如`['A', 'C']`），默认`['A']`
  - `cutoffs`: 碱基个数阈值列表（如`[3, 3]`），默认`[3]`
  - 只有当R1和R2在所有条件下都满足（每个碱基计数都小于对应cutoff）时才保留pair
  - 自动对过滤后的BAM进行排序并建立索引
  - 返回: (原始BAM路径, 过滤BAM路径)

- `parse_index_file(index_file: Path) -> dict`
  - 解析索引文件，返回`{read_id: {base: [indices]}}`字典

- `restore_sequence(read: pysam.AlignedSegment, read_indices: dict) -> bool`
  - 根据索引文件恢复read的序列（保留原始质量值）

- `count_base_in_read(read: pysam.AlignedSegment, target_base: str) -> int`
  - 统计read中指定碱基的个数（考虑链向）

### pileup_PE.py

**Pileup分析模块**（位于`scripts/`目录）

- `get_genome_wide_pileup_paired(bam_path: Path, fasta_path: Path, outputdir: Path, thread: int, prefix: str, tmp_dir: str | None = None, chunk_size: int = 500000, max_depth: int = 100000, allowed_ref_bases: set = None, primary_read: str = 'R1', logger: Optional[logging.Logger] = None) -> Path`
  - 生成双端配对处理的pileup
  - `primary_read`: 'R1' 或 'R2'，指定主要统计哪个read（DNA流程使用R2，RNA流程使用R1）
  - `allowed_ref_bases`: 允许的参考碱基集合（如`{'A', 'C'}`），默认`{'A'}`
  - 使用多进程处理，每个worker写入临时文件，最后合并
  - 返回: pileup文件路径

### call_site_FDR.py

**位点检测模块**（位于`scripts/`目录）

- `call_site_FDR(pileup_file: Path, output_dir: Path, prefix: str, base_type: str, molecule_type: str, min_pair_cov: int = 10, min_base_cov: int = 2, min_base_rate: float = 0.1, fdr_threshold: float = 0.05, logger: Optional[logging.Logger] = None) -> Tuple[Path, Path]`
  - 从pileup调用修饰位点（FDR校正）
  - `base_type`: 'A' 或 'C'
  - `molecule_type`: 'DNA' 或 'RNA'
    - DNA的A是6mA，C是5mC
    - RNA的A是m6A，C是m5C
  - 使用二项检验和Benjamini-Hochberg FDR校正
  - 返回: (位点文件路径, 转化率文件路径)

### Spike_in_analysis.py

**Spike-in分析模块**（位于根目录）

- `analyze_spikein_methylation(spikein_fa: Path, index_file: Path, r1_fastq: Path, r2_fastq: Path, tooldir: Path, outputdir: Path, prefix: str, type: str, threads: int = 8, logger: Optional[logging.Logger] = None) -> Path`
  - 分析spike-in序列的甲基化水平
  - 自动转换spike-in FASTA（根据type）
  - 构建bowtie2索引
  - 转换FASTQ文件
  - 执行mapping（自动排序和索引）
  - 生成pileup
  - 统计目标位点的甲基化水平
  - 返回: 甲基化水平统计结果文件路径

- `parse_index_file(index_file: Path) -> Dict[str, List[int]]`
  - 解析spike-in index文件（格式：序列ID后跟位置列表，用`_`分隔）

- `analyze_pileup_for_positions(pileup_file: Path, index_dict: Dict[str, List[int]], spikein_fa: Path, result_file: Path, logger: Optional[logging.Logger] = None) -> Dict[str, Dict[int, Dict[str, float]]]`
  - 从pileup文件中统计目标位点的甲基化水平

## 命名规范

1. **函数名**：使用小写字母和下划线（snake_case）
2. **参数名**：使用小写字母和下划线（snake_case）
3. **类型参数**：使用大写字母开头的驼峰命名（PascalCase）
4. **常量**：使用全大写字母和下划线（ALL_CAPS）

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

对于根目录的模块（如`pre_process.py`、`download.py`），导入scripts目录的模块：

```python
try:
    from .scripts.run_cmd_and_logging import run_cmd, get_logger
except ImportError:
    from scripts.run_cmd_and_logging import run_cmd, get_logger
```

## 错误处理

所有函数应该：
1. 验证输入参数
2. 提供清晰的错误消息
3. 记录关键操作到日志
4. 验证输出文件是否成功创建

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
