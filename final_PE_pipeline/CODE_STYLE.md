# 代码风格规范

本文档定义了`final_PE_pipeline`项目的统一代码风格规范。

## 文件头格式

所有Python脚本文件应包含以下文件头：

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
模块名称

模块功能简要描述（1-2行）
"""
```

## 导入顺序

导入语句应按以下顺序组织，每组之间用空行分隔：

1. **标准库导入**（按字母顺序）
2. **第三方库导入**（按字母顺序）
3. **本地模块导入**（相对导入优先）

```python
# 标准库
import argparse
import logging
import sys
from pathlib import Path
from typing import Optional, Tuple

# 第三方库
import pysam
from Bio import SeqIO
from Bio.Seq import Seq

# 本地模块
from .run_cmd_and_logging import run_cmd, get_logger
```

## 函数文档字符串

所有公开函数必须包含详细的中文文档字符串，格式如下：

```python
def function_name(
    param1: Type,
    param2: Type = default_value,
    logger: Optional[logging.Logger] = None
) -> ReturnType:
    """
    函数功能简要描述（1行）
    
    详细描述（可选，多行）
    
    参数:
        param1: 参数1描述
        param2: 参数2描述（默认值：xxx）
        logger: logger对象（可选），如果不提供则使用默认logger
    
    返回:
        返回值描述
    
    异常:
        ValueError: 当参数无效时抛出
        FileNotFoundError: 当文件不存在时抛出
    """
```

## 类型提示

- 所有函数参数和返回值都应包含类型提示
- 使用`Optional[Type]`表示可选参数
- 使用`Tuple[Type1, Type2]`表示元组返回类型
- 使用`List[Type]`、`Dict[KeyType, ValueType]`等

## Logger处理

所有函数都应支持`logger`参数：

```python
from .run_cmd_and_logging import get_logger

def my_function(param: str, logger: Optional[logging.Logger] = None):
    """
    函数描述
    """
    if logger is None:
        logger = get_logger()
    
    logger.info("执行操作...")
```

## 命令执行

所有Linux命令必须通过`run_cmd()`执行：

```python
from .run_cmd_and_logging import run_cmd

# 正确
run_cmd("samtools sort -o output.bam input.bam", logger=logger)

# 错误：直接使用subprocess
# subprocess.run("samtools sort ...")  # 不允许
```

## 命名规范

1. **函数名**：使用小写字母和下划线（snake_case）
   - 例如：`extract_longest_isoform`, `build_star_index`

2. **变量名**：使用小写字母和下划线（snake_case）
   - 例如：`genome_fa`, `output_dir`, `r1_fastq`

3. **常量**：使用大写字母和下划线（UPPER_SNAKE_CASE）
   - 例如：`_COMP_BASE_MAP`, `MAX_DEPTH`

4. **类名**：使用大写字母开头的驼峰命名（PascalCase）
   - 例如：`LoggerWriter`

5. **私有函数/变量**：使用单下划线前缀
   - 例如：`_parse_convert_pairs`, `_open_read_handle`

## 注释规范

1. **行内注释**：使用`#`，与代码之间至少一个空格
   ```python
   # 验证输入参数
   if not input_file.exists():
       raise FileNotFoundError(...)
   ```

2. **块注释**：用于解释复杂逻辑
   ```python
   # 处理逻辑：
   # 1. 首先验证输入
   # 2. 然后执行转换
   # 3. 最后验证输出
   ```

3. **文档字符串**：用于函数和类的说明（见上文）

## 错误处理

1. 验证输入参数，提供清晰的错误消息
2. 使用适当的异常类型
3. 记录关键操作到日志

```python
if not input_file.exists():
    raise FileNotFoundError(f"Input file not found: {input_file}")

if threads < 1:
    raise ValueError(f"threads must be >= 1, got {threads}")
```

## 代码组织

1. **函数顺序**：
   - 私有辅助函数在前
   - 公开函数在后
   - 主函数（如果有）在最后

2. **空行使用**：
   - 函数定义之间用2个空行分隔
   - 类定义之间用2个空行分隔
   - 逻辑块之间用1个空行分隔

3. **行长度**：建议不超过100字符，必要时使用括号换行

## 示例

完整的函数示例：

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
示例模块

展示标准代码风格
"""

# 标准库
import logging
from pathlib import Path
from typing import Optional

# 第三方库
import pysam

# 本地模块
from .run_cmd_and_logging import run_cmd, get_logger


def process_bam_file(
    input_bam: Path,
    output_bam: Path,
    threads: int = 8,
    logger: Optional[logging.Logger] = None
) -> Path:
    """
    处理BAM文件
    
    对输入的BAM文件进行排序和索引。
    
    参数:
        input_bam: 输入BAM文件路径
        output_bam: 输出BAM文件路径
        threads: 线程数（默认：8）
        logger: logger对象（可选）
    
    返回:
        输出BAM文件路径
    
    异常:
        FileNotFoundError: 当输入文件不存在时抛出
        ValueError: 当threads < 1时抛出
    """
    if logger is None:
        logger = get_logger()
    
    # 验证输入
    if not input_bam.exists():
        raise FileNotFoundError(f"Input BAM not found: {input_bam}")
    
    if threads < 1:
        raise ValueError(f"threads must be >= 1, got {threads}")
    
    # 确保输出目录存在
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    
    # 执行排序
    logger.info(f"Sorting BAM: {input_bam}")
    cmd = f"samtools sort -@ {threads} -o {output_bam} {input_bam}"
    run_cmd(cmd, logger=logger)
    
    # 验证输出
    if not output_bam.exists():
        raise FileNotFoundError(f"Output BAM was not created: {output_bam}")
    
    logger.info(f"BAM processing completed: {output_bam}")
    return output_bam
```

