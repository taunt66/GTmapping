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


def sort_and_index_bam(
    bam_path: Path,
    logger: Optional[logging.Logger] = None,
    already_sorted: bool = False,
    threads: Optional[int] = None,
) -> Path:
    """
    使用 samtools 对 BAM 文件进行排序并建立索引。

    参数:
        bam_path: 输入 BAM 文件路径（会在原地被排序覆盖）
        logger: logger 对象（可选）
        already_sorted: 如果为 True，则跳过 sort 只建立索引（默认 False）
        threads: samtools sort 使用的线程数（可选，不指定则使用默认线程数）

    返回:
        排序并建立索引后的 BAM 文件路径（与输入路径相同）
    """
    if logger is None:
        logger = logging.getLogger("pipeline")

    bam_path = Path(bam_path)

    # 排序（如有需要）
    if not already_sorted:
        sorted_bam = bam_path.with_suffix(".sorted.bam")
        thread_opt = f"-@ {threads} " if threads else ""
        sort_cmd = f"samtools sort {thread_opt}-o {sorted_bam} {bam_path}"
        logger.info(f"[BAM] Sorting BAM: {bam_path} -> {sorted_bam}")
        run_cmd(sort_cmd, logger=logger)

        if sorted_bam.exists():
            # 删除未排序的原始文件，用排序后的文件覆盖
            try:
                bam_path.unlink()
            except Exception as e:
                logger.warning(f"[BAM] Failed to remove unsorted BAM {bam_path}: {e}")
            sorted_bam.rename(bam_path)
            logger.info(f"[BAM] BAM sorted (in-place): {bam_path}")
        else:
            logger.warning(f"[BAM] Sorted BAM not found after samtools sort: {sorted_bam}")

    # 建立索引
    index_cmd = f"samtools index {bam_path}"
    logger.info(f"[BAM] Indexing BAM: {bam_path}")
    run_cmd(index_cmd, logger=logger)

    index_file = bam_path.with_suffix(bam_path.suffix + ".bai")
    if index_file.exists():
        logger.info(f"[BAM] BAM index file created: {index_file}")
    else:
        logger.warning(f"[BAM] BAM index file not found: {index_file}")

    return bam_path


# 便捷函数：获取默认 logger
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

