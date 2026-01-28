# DNA和RNA流程逻辑说明（默认参数）

本文档详细说明在默认参数下，DNA和RNA两种类型的完整处理流程。

## 默认参数设置

### 通用默认参数
- `min_pair_cov`: 10（位点调用最小pair覆盖度）
- `min_base_cov`: 2（位点调用最小碱基覆盖度）
- `min_base_rate`: 0.1（位点调用最小碱基比率）
- `fdr_threshold`: 0.05（FDR阈值）

### DNA默认参数
- `filter_target_bases`: `['T', 'G']`（过滤目标碱基）
- `filter_cutoffs`: `[3, 3]`（过滤阈值，T和G都使用cutoff=3）

### RNA默认参数
- `filter_target_bases`: `['A', 'C']`（过滤目标碱基）
- `filter_cutoffs`: `[3, 3]`（过滤阈值，A和C都使用cutoff=3）

---

## Step 0: FASTQ转换和索引生成（根据type不同）

**输入：** 原始R1和R2 FASTQ文件

### DNA流程

**转换规则：**
- **R1**: `T->C`, `G->A`
- **R2**: `A->G`, `C->T`

**输出：**
- `{prefix}_r1_converted.fq.gz`: 转换后的R1 FASTQ
- `{prefix}_r2_converted.fq.gz`: 转换后的R2 FASTQ
- `{prefix}_r1_index_T-C_G-A.txt`: R1索引文件（记录原始T和G的位置）
- `{prefix}_r2_index_A-G_C-T.txt`: R2索引文件（记录原始A和C的位置）

### RNA流程

**转换规则：**
- **R1**: `A->G`, `C->T`
- **R2**: `T->C`, `G->A`

**输出：**
- `{prefix}_r1_converted.fq.gz`: 转换后的R1 FASTQ
- `{prefix}_r2_converted.fq.gz`: 转换后的R2 FASTQ
- `{prefix}_r1_index_A-G_C-T.txt`: R1索引文件（记录原始A和C的位置）
- `{prefix}_r2_index_T-C_G-A.txt`: R2索引文件（记录原始T和G的位置）

**说明：** 转换后的FASTQ用于后续mapping，索引文件用于后续还原原始BAM。

---

## Step 1: 文件准备（根据type不同）

### DNA流程

**需要准备的文件：**
1. **基因组TC_GA转换文件** (`{genome}_TC_GAconvert.fa`)
   - 转换规则：`T->C`, `G->A`
   - 用途：第一次mapping的参考基因组

2. **基因组AG_CT转换文件** (`{genome}_AG_CTconvert.fa`)
   - 转换规则：`A->G`, `C->T`
   - 用途：第二次mapping的参考基因组

3. **STAR索引**
   - `{genome}_TC_GA_STAR_index`: TC_GA基因组的STAR索引
   - `{genome}_AG_CT_STAR_index`: AG_CT基因组的STAR索引
   - 注意：DNA流程不使用GTF文件构建索引（no-splice alignment）

### RNA流程

**需要准备的文件：**
1. **最长转录本文件**
   - `{genome}_longest_transcript.fa`: 最长转录本FASTA
   - `{genome}_longest.trans`: 最长转录本信息TSV（包含exon信息）

2. **基因组AG_CT转换文件** (`{genome}_AG_CTconvert.fa`)
   - 转换规则：`A->G`, `C->T`
   - 用途：第一次mapping和转录组转换的参考基因组

3. **基因组TC_GA转换文件** (`{genome}_TC_GAconvert.fa`)
   - 转换规则：`T->C`, `G->A`
   - 用途：第二次mapping的参考基因组

4. **最长转录本AG_CT转换文件** (`{genome}_longest_transcript_AG_CTconvert.fa`)
   - 转换规则：`A->G`, `C->T`
   - 用途：第三次mapping的参考转录组

5. **索引**
   - `{genome}_AG_CT_STAR_index`: AG_CT基因组的STAR索引（不使用GTF）
   - `{genome}_TC_GA_STAR_index`: TC_GA基因组的STAR索引（不使用GTF）
   - `{genome}_longest_transcript_AG_CT_bowtie2_index`: 最长转录本AG_CT的bowtie2索引

---

## Step 2: Mapping流程（根据type不同）

### DNA流程（两步mapping）

#### Step 2a: 第一次mapping - STAR基因组TC_GA（R1正链）

**输入：** 转换后的R1和R2 FASTQ（Step 0的输出）

**Mapping参数：**
- 工具：STAR
- 参考基因组：`{genome}_TC_GAconvert.fa`（T->C, G->A转换）
- 类型：DNA（no-splice alignment）
- 参数：`--alignIntronMax 0`（禁用splicing）

**过滤条件：**
- 提取R1映射到**正链**的unique primary pairs

**输出：**
- `{prefix}.star_genome_TC_GA.filtered.bam`: 过滤后的BAM（R1正链）
- `{prefix}.remaining_1.r1.fq.gz`: 未选中的R1 reads
- `{prefix}.remaining_1.r2.fq.gz`: 未选中的R2 reads

#### Step 2b: 第二次mapping - STAR基因组AG_CT（R1负链）

**输入：** Step 2a剩余的R1和R2 FASTQ

**Mapping参数：**
- 工具：STAR
- 参考基因组：`{genome}_AG_CTconvert.fa`（A->G, C->T转换）
- 类型：DNA（no-splice alignment）

**过滤条件：**
- 提取R1映射到**负链**的unique primary pairs

**输出：**
- `{prefix}.star_genome_AG_CT.filtered.bam`: 过滤后的BAM（R1负链）

#### Step 2c: 合并BAM文件

**输入：** Step 2a和Step 2b的两个过滤BAM

**输出：**
- `{prefix}.merged.sorted.bam`: 合并并排序后的BAM

---

### RNA流程（三步mapping）

#### Step 2a: 第一次mapping - STAR基因组AG_CT（R1正链）

**输入：** 转换后的R1和R2 FASTQ（Step 0的输出）

**Mapping参数：**
- 工具：STAR
- 参考基因组：`{genome}_AG_CTconvert.fa`（A->G, C->T转换）
- 类型：RNA（splicing-aware alignment）

**过滤条件：**
- 提取R1映射到**正链**的unique primary pairs

**输出：**
- `{prefix}.star_genome_AG_CT.filtered.bam`: 过滤后的BAM（R1正链）
- `{prefix}.remaining_1.r1.fq.gz`: 未选中的R1 reads
- `{prefix}.remaining_1.r2.fq.gz`: 未选中的R2 reads

#### Step 2b: 第二次mapping - STAR基因组TC_GA（R1负链）

**输入：** Step 2a剩余的R1和R2 FASTQ

**Mapping参数：**
- 工具：STAR
- 参考基因组：`{genome}_TC_GAconvert.fa`（T->C, G->A转换）
- 类型：RNA（splicing-aware alignment）

**过滤条件：**
- 提取R1映射到**负链**的unique primary pairs

**输出：**
- `{prefix}.star_genome_TC_GA.filtered.bam`: 过滤后的BAM（R1负链）
- `{prefix}.remaining_2.r1.fq.gz`: 未选中的R1 reads
- `{prefix}.remaining_2.r2.fq.gz`: 未选中的R2 reads

#### Step 2c: 第三次mapping - bowtie2转录组AG_CT（R1正链）

**输入：** Step 2b剩余的R1和R2 FASTQ

**Mapping参数：**
- 工具：bowtie2
- 参考转录组：`{genome}_longest_transcript_AG_CTconvert.fa`（A->G, C->T转换）

**过滤条件：**
- 提取R1映射到**正链**的unique primary pairs

**输出：**
- `{prefix}.bowtie2_transcriptome_AG_CT.filtered.bam`: 过滤后的BAM（转录组坐标）

#### Step 2d: 转录组到基因组坐标转换

**输入：** Step 2c的转录组BAM

**转换参数：**
- 参考基因组：`{genome}_AG_CTconvert.fa`
- 转录本信息：`{genome}_longest.trans`

**输出：**
- `{prefix}.bowtie2_transcriptome_AG_CT.genome.bam`: 基因组坐标的BAM

#### Step 2e: 合并三个BAM文件

**输入：** Step 2a、Step 2b和Step 2d的三个BAM

**输出：**
- `{prefix}.merged.sorted.bam`: 合并并排序后的BAM

---

## Step 8: 还原原始BAM并过滤（根据type不同）

**输入：** 合并后的BAM（`{prefix}.merged.sorted.bam`）

**还原过程：**
- 使用Step 0生成的索引文件（`r1_index_file`和`r2_index_file`）
- 将BAM中的转换碱基还原为原始碱基

**过滤条件（默认参数）：**

### DNA流程
- `filter_target_bases`: `['T', 'G']`
- `filter_cutoffs`: `[3, 3]`
- **过滤逻辑：** 只有当R1和R2的T计数都 < 3 **且** G计数都 < 3时，才保留该pair

### RNA流程
- `filter_target_bases`: `['A', 'C']`
- `filter_cutoffs`: `[3, 3]`
- **过滤逻辑：** 只有当R1和R2的A计数都 < 3 **且** C计数都 < 3时，才保留该pair

**输出：**
- `{prefix}.origin.bam`: 还原后的原始BAM
- `{prefix}.filtered.bam`: 过滤后的BAM（用于后续pileup）

---

## Step 9: Pileup生成（根据type不同）

**输入：** 过滤后的BAM（`{prefix}.filtered.bam`）

**参考基因组：** 原始基因组FASTA（**未转换**的`genome_fa`）

**Pileup参数（根据type不同）：**

### DNA流程
- `primary_read`: `'R2'`（主要使用R2进行计数）
- `allowed_ref_bases`: `{'A', 'C'}`（统计A和C的pileup信息）

### RNA流程
- `primary_read`: `'R1'`（主要使用R1进行计数）
- `allowed_ref_bases`: `{'A', 'C'}`（统计A和C的pileup信息）

**输出：**
- `{prefix}.pileup`: Pileup结果文件

---

## Step 10: 修饰位点调用（根据type不同）

**输入：** Pileup文件（`{prefix}.pileup`）

**调用参数（默认）：**
- `min_pair_cov`: 10
- `min_base_cov`: 2
- `min_base_rate`: 0.1
- `fdr_threshold`: 0.05

### DNA流程

**调用的修饰类型：**
1. **6mA**（A位点）
   - `base_type`: `'A'`
   - `molecule_type`: `'DNA'`

2. **5mC**（C位点）
   - `base_type`: `'C'`
   - `molecule_type`: `'DNA'`

**输出：**
- `{prefix}.6mA.sites.tsv`: 6mA修饰位点
- `{prefix}.5mC.sites.tsv`: 5mC修饰位点
- `{prefix}.6mA.conversion_rate.tsv`: 6mA转换率
- `{prefix}.5mC.conversion_rate.tsv`: 5mC转换率

### RNA流程

**调用的修饰类型：**
1. **m6A**（A位点）
   - `base_type`: `'A'`
   - `molecule_type`: `'RNA'`

2. **m5C**（C位点）
   - `base_type`: `'C'`
   - `molecule_type`: `'RNA'`

**输出：**
- `{prefix}.m6A.sites.tsv`: m6A修饰位点
- `{prefix}.m5C.sites.tsv`: m5C修饰位点
- `{prefix}.m6A.conversion_rate.tsv`: m6A转换率
- `{prefix}.m5C.conversion_rate.tsv`: m5C转换率

---

## 总结对比表

| 步骤 | DNA流程 | RNA流程 |
|------|---------|--------|
| **Step 0: FASTQ转换** | R1: T->C, G->A<br>R2: A->G, C->T | R1: A->G, C->T<br>R2: T->C, G->A |
| **Step 1: 文件准备** | 2个基因组转换文件<br>2个STAR索引 | 最长转录本+3个转换文件<br>3个索引（2个STAR+1个bowtie2） |
| **Step 2: Mapping** | 2步：TC_GA正链 + AG_CT负链 | 3步：AG_CT正链 + TC_GA负链 + 转录组正链 |
| **Step 2合并** | 合并2个BAM | 合并3个BAM（含转录组转换） |
| **Step 8: 过滤** | 过滤T和G（cutoff=3） | 过滤A和C（cutoff=3） |
| **Step 9: Pileup** | primary_read='R2'<br>统计A和C | primary_read='R1'<br>统计A和C |
| **Step 10: 位点调用** | 6mA和5mC | m6A和m5C |

---

## 关键设计原理

1. **FASTQ转换的目的：** 通过碱基转换，将修饰位点（如m6A的A->G）转换为可检测的突变，便于后续mapping和检测。

2. **多步mapping的目的：** 
   - DNA：分别捕获正链和负链的reads
   - RNA：分别捕获正链、负链和转录组reads，提高覆盖度

3. **过滤的目的：** 去除含有过多转换碱基的reads，这些可能是假阳性或低质量reads。

4. **Pileup使用原始基因组的原因：** Pileup需要统计原始碱基的覆盖情况，而不是转换后的碱基。

5. **DNA和RNA的差异：**
   - DNA使用no-splice alignment（`--alignIntronMax 0`）
   - RNA使用splicing-aware alignment
   - DNA按R2计数，RNA按R1计数（可能与链向和修饰检测策略有关）

