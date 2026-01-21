"""
final_PE_pipeline.scripts

双端测序数据处理流程的功能模块

所有模块都提供函数接口，支持统一的参数风格和日志记录。
"""

# 版本信息
__version__ = "2.0.0"

# 导入主要功能（可选，用于方便访问）
from .run_cmd_and_logging import setup_project_logger, get_logger, run_cmd

__all__ = [
    "setup_project_logger",
    "get_logger", 
    "run_cmd",
]

