# extract_target_aln_dump_others.py 生成文件说明

## 文件生成流程

`extract_target_aln_and_dump_others` 函数在提取符合条件的比对结果时，会生成以下文件：

### 第一步：从BAM提取未选中的reads（临时文件）

**文件1: `{prefix}.remaining_1.r1.fq.gz`**
- **含义**：从BAM文件中提取的**未选中**的R1 reads的原始FASTQ文件
- **内容**：包括unmapped、链向不符合要求、多重比对的R1 reads
- **特点**：
  - 这些reads**未经过配对和排序**
  - read的顺序可能与R2文件不一致
  - 可能包含单端reads（mate未找到的情况）
- **用途**：作为`seqkit pair`的输入文件之一

**文件2: `{prefix}.remaining_1.r2.fq.gz`**
- **含义**：从BAM文件中提取的**未选中**的R2 reads的原始FASTQ文件
- **内容**：包括unmapped、链向不符合要求、多重比对的R2 reads
- **特点**：
  - 这些reads**未经过配对和排序**
  - read的顺序可能与R1文件不一致
  - 可能包含单端reads（mate未找到的情况）
- **用途**：作为`seqkit pair`的输入文件之一

### 第二步：使用seqkit pair重新配对和排序

**文件3: `{prefix}.remaining_1.r1.paired.fq.gz`**
- **含义**：经过`seqkit pair`重新配对和排序后的R1 FASTQ文件
- **内容**：与R2文件**配对成功**的R1 reads
- **特点**：
  - **已配对**：每个R1都有对应的R2
  - **已排序**：R1和R2的read顺序完全一致（按read name排序）
  - **已过滤**：只包含能够配对的reads，单端reads已被过滤
- **用途**：用于后续的mapping步骤（作为下一次mapping的输入）

**文件4: `{prefix}.remaining_1.r2.paired.fq.gz`**
- **含义**：经过`seqkit pair`重新配对和排序后的R2 FASTQ文件
- **内容**：与R1文件**配对成功**的R2 reads
- **特点**：
  - **已配对**：每个R2都有对应的R1
  - **已排序**：R1和R2的read顺序完全一致（按read name排序）
  - **已过滤**：只包含能够配对的reads，单端reads已被过滤
- **用途**：用于后续的mapping步骤（作为下一次mapping的输入）

## 文件大小差异说明

从您提供的文件大小来看：
- `H_L1.remaining_1.r1.fq.gz`: 10M（未配对的原始文件）
- `H_L1.remaining_1.r1.paired.fq.gz`: 12M（已配对的文件）
- `H_L1.remaining_1.r2.fq.gz`: 9.7M（未配对的原始文件）
- `H_L1.remaining_1.r2.paired.fq.gz`: 12M（已配对的文件）

**为什么`.paired`文件比原始文件大？**

这是因为：
1. **原始文件**（`.r1.fq.gz`和`.r2.fq.gz`）包含：
   - 所有未选中的reads（包括单端reads）
   - 未排序，可能有重复或缺失的reads

2. **配对文件**（`.paired.fq.gz`）包含：
   - 只包含能够配对的reads（R1和R2都存在）
   - 已排序，确保R1和R2的顺序一致
   - `seqkit pair`可能会重新组织数据，导致压缩率不同

**注意**：文件大小的差异是正常的，因为：
- 压缩率可能不同（排序后的数据压缩率可能更高或更低）
- 配对过程会过滤掉单端reads
- 数据重新组织可能影响gzip压缩效果

## 工作流程总结

```
输入BAM文件
    ↓
提取符合条件的pairs → 输出到 {prefix}.filtered.bam
    ↓
提取未选中的reads → 输出到 {prefix}.remaining_1.r1.fq.gz 和 .r2.fq.gz（未配对）
    ↓
seqkit pair 重新配对和排序
    ↓
输出配对后的reads → {prefix}.remaining_1.r1.paired.fq.gz 和 .r2.paired.fq.gz（已配对）
    ↓
返回：filtered.bam, r1.paired.fq.gz, r2.paired.fq.gz
```

## 返回值说明

函数返回三个文件：
1. `output_bam`：符合条件的比对结果BAM文件
2. `repaired_r1`：配对后的R1 FASTQ文件（`.paired.fq.gz`）
3. `repaired_r2`：配对后的R2 FASTQ文件（`.paired.fq.gz`）

**注意**：原始的`.r1.fq.gz`和`.r2.fq.gz`文件（未配对的）会在生成配对文件后**自动删除**，只保留`.paired.fq.gz`文件。

