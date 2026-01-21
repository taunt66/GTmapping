#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_cmd_and_logging.py

提供项目级别的日志记录功能：
1. 设置全局logger，捕获Python的stderr和stdout
2. 提供运行shell命令并记录输出的函数
"""

import logging
import sys
import subprocess
from pathlib import Path
from typing import Optional, List
from datetime import datetime
import io


class LoggerWriter:
    """将stdout/stderr重定向到logger的类"""
    def __init__(self, logger, level):
        self.logger = logger
        self.level = level
        self.buffer = io.StringIO()
    
    def write(self, message):
        if message.strip():  # 忽略空行
            self.buffer.write(message)
            # 检查是否有换行符，如果有则立即刷新
            if '\n' in message:
                self.flush()
    
    def flush(self):
        if self.buffer.getvalue():
            self.logger.log(self.level, self.buffer.getvalue().rstrip())
            self.buffer = io.StringIO()


def setup_project_logger(log_dir: Path, log_name: str = "pipeline.log",
                        level: int = logging.INFO) -> logging.Logger:
    """
    设置项目级别的logger，捕获Python的stderr和stdout
    
    参数:
        log_dir: 日志文件目录
        log_name: 日志文件名（默认：pipeline.log）
        level: 日志级别（默认：INFO）
    
    返回:
        logger对象
    """
    # 确保日志目录存在
    log_dir.mkdir(parents=True, exist_ok=True)
    
    # 创建日志文件路径（带时间戳）
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"{log_name}.{timestamp}"
    
    # 配置logger
    logger = logging.getLogger("pipeline")
    logger.setLevel(level)
    
    # 避免重复添加handler
    if logger.handlers:
        return logger
    
    # 文件handler（记录所有日志）
    file_handler = logging.FileHandler(log_file, mode='a', encoding='utf-8')
    file_handler.setLevel(level)
    file_format = logging.Formatter(
        '%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    file_handler.setFormatter(file_format)
    logger.addHandler(file_handler)
    
    # 控制台handler（只显示INFO及以上级别）
    console_handler = logging.StreamHandler(sys.__stdout__)
    console_handler.setLevel(logging.INFO)
    console_format = logging.Formatter(
        '[%(levelname)s] %(message)s'
    )
    console_handler.setFormatter(console_format)
    logger.addHandler(console_handler)
    
    # 重定向stdout和stderr到logger
    sys.stdout = LoggerWriter(logger, logging.INFO)
    sys.stderr = LoggerWriter(logger, logging.ERROR)
    
    logger.info(f"Logger initialized. Log file: {log_file}")
    logger.info("=" * 80)
    
    return logger


def run_cmd(cmd: str, logger: Optional[logging.Logger] = None,
            shell: bool = True, check: bool = True,
            cwd: Optional[Path] = None, env: Optional[dict] = None) -> subprocess.CompletedProcess:
    """
    运行shell命令并将stderr和stdout重定向到logger
    
    参数:
        cmd: 要执行的命令（字符串）
        logger: logger对象，如果为None则使用默认logger
        shell: 是否使用shell执行（默认：True）
        check: 如果命令失败是否抛出异常（默认：True）
        cwd: 工作目录（可选）
        env: 环境变量字典（可选）
    
    返回:
        subprocess.CompletedProcess对象
    
    示例:
        run_cmd("ls -l", logger)
        run_cmd("samtools view input.bam", logger, check=False)
    """
    if logger is None:
        logger = logging.getLogger("pipeline")
    
    # 记录命令
    logger.info(f"[CMD] {cmd}")
    if cwd:
        logger.info(f"[CWD] {cwd}")
    
    try:
        # 执行命令
        result = subprocess.run(
            cmd,
            shell=shell,
            check=check,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,  # 将stderr合并到stdout
            text=True,
            cwd=str(cwd) if cwd else None,
            env=env
        )
        
        # 记录输出
        if result.stdout:
            for line in result.stdout.split('\n'):
                if line.strip():
                    logger.info(f"[OUT] {line}")
        
        # 记录返回码
        if result.returncode == 0:
            logger.info(f"[SUCCESS] Command completed with return code {result.returncode}")
        else:
            logger.warning(f"[WARNING] Command completed with return code {result.returncode}")
        
        return result
        
    except subprocess.CalledProcessError as e:
        logger.error(f"[ERROR] Command failed with return code {e.returncode}")
        if e.stdout:
            for line in e.stdout.split('\n'):
                if line.strip():
                    logger.error(f"[OUT] {line}")
        raise
    except Exception as e:
        logger.error(f"[ERROR] Command execution failed: {e}")
        raise


def run_cmd_list(cmd_list: List[str], logger: Optional[logging.Logger] = None,
                 shell: bool = True, check: bool = True,
                 cwd: Optional[Path] = None, env: Optional[dict] = None) -> List[subprocess.CompletedProcess]:
    """
    运行多个shell命令（顺序执行）
    
    参数:
        cmd_list: 命令列表
        logger: logger对象
        shell: 是否使用shell执行
        check: 如果命令失败是否抛出异常
        cwd: 工作目录（可选）
        env: 环境变量字典（可选）
    
    返回:
        结果列表
    """
    if logger is None:
        logger = logging.getLogger("pipeline")
    
    results = []
    for i, cmd in enumerate(cmd_list, 1):
        logger.info(f"[BATCH] Running command {i}/{len(cmd_list)}")
        result = run_cmd(cmd, logger, shell, check, cwd, env)
        results.append(result)
    
    return results


# 便捷函数：获取默认logger
def get_logger() -> logging.Logger:
    """获取默认的pipeline logger"""
    return logging.getLogger("pipeline")


if __name__ == "__main__":
    # 示例用法
    from pathlib import Path
    
    # 设置logger
    log_dir = Path("./logs")
    logger = setup_project_logger(log_dir, "test_pipeline.log")
    
    # 测试Python输出捕获
    print("This will be captured by logger")
    sys.stderr.write("This error will also be captured\n")
    
    # 测试命令执行
    run_cmd("echo 'Hello from shell command'", logger)
    run_cmd("ls -lh", logger, check=False)
    
    logger.info("=" * 80)
    logger.info("Test completed")

