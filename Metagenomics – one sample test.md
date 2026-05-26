# Metagenomics – one sample test

**项目名称**: Livestock_farm_1 Metagenomics 
**操作日期**: 2025-12
**操作人员**: Shengjie Xu

[TOC]

## 1. Data Verification——md5sum

为了确保下载数据的完整性，我们需要计算本地文件的 MD5 值，并与数据库中提供的md5sum.txt进行比对

```bash
# 1. 进入你的数据目录
cd /ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/1_rawdata
# 2. 批量计算所有 .gz 文件的 MD5 值，并保存到文件
find . -type f -name "*.gz" -print0 | xargs -0 md5sum > md5_calculated.txt
```

创建 `check_txt_vs_txt.py` 脚本，用于批量比对

```bash
nano check_txt_vs_txt.py
```

比对md5_calculated.txt与md5sum.txt

```python
import os

# --- 配置文件名 ---
official_file = "md5sum.txt"       # 官方文件 (Standard)
my_file = "md5_calculated.txt"     # 本地计算文件 (Calculated)

def load_md5_file(filepath):
    data = {}
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                md5 = parts[0]
                # 只保留文件名，忽略前面的路径
                fname = os.path.basename(parts[1])
                data[fname] = md5
    return data

print(f"Loading {official_file} ...")
official_dict = load_md5_file(official_file)

print(f"Loading {my_file} ...")
my_dict = load_md5_file(my_file)

# --- 开始比对 ---
match_count = 0
mismatch = []
missing = []

for fname, true_md5 in official_dict.items():
    if fname in my_dict:
        if my_dict[fname] == true_md5:
            match_count += 1
        else:
            mismatch.append(fname)
    else:
        missing.append(fname)

print("-" * 30)
print(f"✅ 校验通过: {match_count} 个文件")

if mismatch:
    print(f"❌ MD5 不匹配: {len(mismatch)} 个文件")
    for f in mismatch: print(f"   - {f}")
elif missing:
    print(f"⚠️  文件缺失: {len(missing)} 个文件")
else:
    print("🎉 完美！所有文件校验通过 (All Passed)！")
```

按照提示保存并关闭文本编辑器后，运行`check_txt_vs_txt.py` 脚本

```bash
python check_txt_vs_txt.py
```

若 print("🎉 完美...)则可以进行Step 2 质控环节；否则需要重新检查下载数据，是否因网络中断等问题导致的文件缺失。

---

## 2. QC

### 前期环境准备

- [ ] 已安装fastp——environment location： `/share/home/u24149/miniconda3/envs/fastp`

- [ ] 已安装bbmap——environment location：`/share/home/u24149/miniconda3/envs/bbmap

- [ ] 在/share/home/u24149/miniconda3/envs/bbmap/opt/bbmap-39.52-0/resources目录下，下载了一系列reference_sequence，激活 bbmap 环境即可使用。

  ![image-20251208102605226](/Users/xu.sj/Library/Application Support/typora-user-images/image-20251208102605226.png)

### Test: 以CRR1071866为例

```bash
cd /ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1 #切换到自己工作目录下
mkdir -p ./2_qc/test #-p：自动创建所有不存在的上级目录，若目录已存在也不会报错
cd ./2_qc/test
mkdir 1_fastqc_raw_data 2_trim_ployG 3_fastqc_trim_polyG 4_trim_data_bbduk 5_fastqc_quality_data #分步质控
```

#### 总体思路：

| 文件夹名称/步骤       | 主要作用                                                     | 输入数据来源                                                 | 逻辑目的                                                     |
| :-------------------- | :----------------------------------------------------------- | :----------------------------------------------------------- | :----------------------------------------------------------- |
| 1_fastqc_raw_data     | 对raw reads进行 FastQC 质量评估；生成 HTML/zip 报告，检查 Q-score、GC 偏倚、adapter 污染等。 | 原始 FASTQ 文件：这里先用CRR1071866测序raw文件作为例子展示。 | 基准 QC：识别原始问题（如低质量尾巴、adapter），不修改数据。只评估。 |
| 2_trim_ployG          | 使用 fastp初步修剪 polyG/polyX 尾巴。去除 3' 端重复序列，并轻度过滤低质量。 | 原始 FASTQ 文件                                              | 针对性修剪：基于阶段1报告，专注 polyG 问题。                 |
| 3_fastqc_trim_polyG   | 对第二步修剪后的数据进行 FastQC 评估。比较前后变化，确认修剪效果。 | 2_trim_ployG 的输出                                          | 验证修剪：检查是否改善了 Q-score 和 overrepresented sequences。迭代优化。 |
| 4_trim_data_bbduk     | 使用 BBduk（BBMap 套件）高级修剪：去除 adapter、phiX 污染、短 reads 等，即进一步过滤。 | 3_fastqc_trim_polyG 的数据                                   | 深度修剪：处理第三步报告中剩余问题（如污染）。BBduk 高效，适合大文件。 |
| 5_fastqc_quality_data | 对最终修剪数据进行 FastQC，作为最终 QC，并生成汇总报告。     | 4_trim_data_bbduk 的输出                                     | 终检：确保数据 ready for downstream（如组装 MEGAHIT、分类 Kraken）。 |

```bash
cd ./1_fastqc_raw_data
nano 1_fastqc_raw_data.sh #或vim，按照个人习惯
```

#### 第一步：Raw_data FastQC

进入`1_fastqc_raw_data.sh`，输入以下script

```bash
#!/bin/bash

INPUT_DIR="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/1_rawdata/CRR1071866" #根据自己文件所在目录进行修改

# 运行 FastQC
# -t 2: 使用2个线程
# -o .: 结果输出到当前目录
# *.fastq*: 自动匹配目录下所有的 fastq 或 fastq.gz 文件

fastqc -t 2 -f fastq -o . ${INPUT_DIR}/*.fq.gz
```

保存并退出文本编辑器，运行脚本

```bash
bash 1_fastqc_raw_data.sh
#随后生成两个ZIP和两个HTML文件（双端）
```

打开本地Terminal，下载并查看.html文件

```bash
scp u24149@logini.tongji.edu.cn:"/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/2_qc/test/1_fastqc_raw_data/*.html" /Users/xu.sj/Desktop #选择本地下载位置
#输入服务器密码
```

![image-20251205140426378](/Users/xu.sj/Library/Application Support/typora-user-images/image-20251205140426378.png)

以一端为例，FastQC 的测序数据质量报告显示：**红色曲线（实际测序序列的 GC 分布）与蓝色曲线（理论正态分布）严重偏离**。环境样本包含多种微生物（细菌、真菌、病毒、古菌），每个物种 GC 含量差异大（细菌 25-75%，真核 40-60%）。导致整体 reads GC 分布宽广、非正态，属于正常现象。

#### 第二步：Trim polyG

```bash
cd ../2_trim_polyG
nano 2_trim_raw_data_polyG.sh
```

进入`2_trim_raw_data_polyG.sh`，输入以下script

```bash
#!/bin/bash
# Fastp trimming script for sample CRR1071866
set -e # Exit if any command fails

# Define input/output
R1_IN="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/1_rawdata/CRR1071866/CRR1071866_f1.fq.gz"
R2_IN="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/1_rawdata/CRR1071866/CRR1071866_r2.fq.gz"
R1_OUT="CRR1071866_f1_trim_polyG.fastq.gz"
R2_OUT="CRR1071866_r2_trim_polyG.fastq.gz"

# 日志文件名加上样本名，方便查看
LOG="fastp_CRR1071866_$(date +%Y%m%d_%H%M).log"

# Activate fastp
source /share/home/u24149/miniconda3/etc/profile.d/conda.sh
conda activate fastp

# Run fastp
echo "[$(date)] Starting fastp trimming for CRR1071866..." | tee "$LOG"

fastp \
  -i "$R1_IN" -I "$R2_IN" \
  -o "$R1_OUT" -O "$R2_OUT" \
  --detect_adapter_for_pe \ #自动检测并修剪 paired-end 数据中的 adapter 序列
  --trim_poly_g \ #修剪 3' 端 polyG 序列
  --trim_poly_x \ #修剪 3' 端任意 polyX 序列（X 为 A/C/G/T，低复杂度区域）
  --cut_tail --cut_mean_quality 20 --cut_window_size 4 \ #启用尾部质量修剪：从 3' 端切除低质量碱基。修剪阈值：窗口内平均 Q-score <20 时切除。滑动窗口大小：4 bp 窗口评估质量。
  --length_required 75 \ #过滤太短序列：最终长度 <75 bp 的 reads 丢弃。
  --thread 8 \
  --html CRR1071866_fastp.html \
  --json CRR1071866_fastp.json 2>&1 | tee -a "$LOG"

echo "[$(date)] fastp finished successfully." | tee -a "$LOG"
```

保存并退出文本编辑器，运行脚本

```bash
bash 2_trim_raw_data_polyG.sh
```

打开本地Terminal，下载并查看.html文件

```bash
scp u24149@logini.tongji.edu.cn:"/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/2_qc/test/2_trim_polyG/*.html" /Users/xu.sj/Desktop #选择本地下载位置
#输入服务器密码
```

![image-20251205154538205](/Users/xu.sj/Library/Application Support/typora-user-images/image-20251205154538205.png)

fastp 报告显示GC 含量稳定 ~63%，低质量/短 reads 极少。

#### 第三步：Trim_raw_data FastQC

```bash
cd ../3_fastqc_trim_polyG
nano 3_fastqc_trim_data.sh
```

进入`3_fastqc_trim_data.sh`，输入以下script

```bash
#!/bin/bash

INPUT_DIR="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/2_qc/test/2_trim_polyG" #根据自己文件所在目录进行修改

fastqc -t 2 -f fastq -o . ${INPUT_DIR}/*.fastq.gz
```

保存并退出文本编辑器，运行脚本

```bash
bash 3_fastqc_trim_data.sh
```

打开本地Terminal，下载并查看.html文件

```bash
scp u24149@logini.tongji.edu.cn:"/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/2_qc/test/3_fastqc_trim_polyG/*.html" /Users/xu.sj/Desktop
#输入服务器密码
```

![image-20251205165547393](/Users/xu.sj/Library/Application Support/typora-user-images/image-20251205165547393.png)

Trimming 前后只损失了 **295 条** Reads，说明原始数据的质量极高，绝大部分 Reads 都保留下来了。由于在 `fastp` 脚本中设置了 --length_required 75，Sequence length 分布变为 **75-150 bp**。

#### 第四步：BBDuk Deep Cleaning

强制修剪末端（ftm=5） → 强制修剪左端10bp（ftl=10） → 去除接头（Adapters） → 过滤PhiX污染 → 质量过滤


```bash
cd ../4_trim_data_bbduk
nano 4_trim_data_bbduk.sh
```

进入`4_trim_data_bbduk.sh`，输入以下script

```bash
#!/bin/bash
# BBDuk Deep Cleaning Pipeline for sample CRR1071866

set -euo pipefail # 遇到错误立即停止

# ================= 1. 环境设置 =================
# 激活 conda 环境
source /share/home/u24149/miniconda3/etc/profile.d/conda.sh
conda activate bbmap

# ================= 2. 输入输出定义 =================
R1_IN="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/2_qc/test/2_trim_polyG/CRR1071866_f1_trim_polyG.fastq.gz"
R2_IN="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/2_qc/test/2_trim_polyG/CRR1071866_r2_trim_polyG.fastq.gz"

# 参考文件路径
ADAPTERS="adapters" 
PHIX="phix"

# Final outputs
R1_OUT="CRR1071866_qc.R1.fastq"
R2_OUT="CRR1071866_qc.R2.fastq"
STATS_TXT="CRR1071866_bbduk_stats.txt"

# 内存设置 
MEM="10g"

# 中间文件定义
TMP1_R1="tmp_1_ftm_R1.fq"
TMP1_R2="tmp_1_ftm_R2.fq"
TMP2_R1="tmp_2_ftl_R1.fq"
TMP2_R2="tmp_2_ftl_R2.fq"
TMP3_R1="tmp_3_ada_R1.fq"
TMP3_R2="tmp_3_ada_R2.fq"
TMP4_R1="tmp_4_phix_R1.fq"
TMP4_R2="tmp_4_phix_R2.fq"
PHIX_MATCHED="CRR1071866_matched_phix.fq.gz"

# 是否清理中间文件
CLEAN_INTERMEDIATES=false

# 日志
LOG="bbduk_CRR1071866_$(date +%Y%m%d_%H%M).log"

# ================= 3. 开始运行 =================
echo "[$(date)] Starting BBDuk pipeline for CRR1071866..." | tee "$LOG"

# 检查输入文件是否存在
if [ ! -f "$R1_IN" ]; then
    echo "Error: Input file $R1_IN not found!" | tee -a "$LOG"
    exit 1
fi

# Step 1: Trim Last Base (ftm=5) - 修剪末端以适配 5 的倍数
echo "--- Step 1: Force Trim Modulo (ftm=5) ---" | tee -a "$LOG"
bbduk.sh -Xmx"$MEM" in="$R1_IN" in2="$R2_IN" out="$TMP1_R1" out2="$TMP1_R2" \
    ftm=5 overwrite=t 2>&1 | tee -a "$LOG"

# Step 2: Trim Left 10bp (ftl=10) - 切除头部 10bp：解决在CRR1071866_f1_trim_polyG.fastq.gz文件中，Per base sequence content一栏，前10bp不稳定的情况（常见于测序起始偏差）
echo "--- Step 2: Force Trim Left 10bp ---" | tee -a "$LOG"
bbduk.sh -Xmx"$MEM" in="$TMP1_R1" in2="$TMP1_R2" out="$TMP2_R1" out2="$TMP2_R2" \
    ftl=10 overwrite=t 2>&1 | tee -a "$LOG"

# Step 3: Adapter Trim - 去接头
echo "--- Step 3: Adapter Trimming ---" | tee -a "$LOG"
bbduk.sh -Xmx"$MEM" in="$TMP2_R1" in2="$TMP2_R2" out="$TMP3_R1" out2="$TMP3_R2" \
    ref="$ADAPTERS" tbo tpe k=23 mink=11 hdist=1 ktrim=r overwrite=t 2>&1 | tee -a "$LOG"

# Step 4: PhiX Filter - 去除 PhiX 污染（常见于NGS库构建中的spike-in对照）
echo "--- Step 4: PhiX Filtering ---" | tee -a "$LOG"
bbduk.sh -Xmx"$MEM" in="$TMP3_R1" in2="$TMP3_R2" out="$TMP4_R1" out2="$TMP4_R2" \
    outm="$PHIX_MATCHED" ref="$PHIX" k=31 hdist=1 stats="$STATS_TXT" overwrite=t 2>&1 | tee -a "$LOG"

# Step 5: Quality Filter - 质量过滤并压缩输出
echo "--- Step 5: Quality Filtering & Compression ---" | tee -a "$LOG"
bbduk.sh -Xmx"$MEM" in="$TMP4_R1" in2="$TMP4_R2" out="$R1_OUT" out2="$R2_OUT" \
    qtrim=rl trimq=15 minlength=30 overwrite=t 2>&1 | tee -a "$LOG"

# ================= 4. 清理与完成 =================

if [ "$CLEAN_INTERMEDIATES" = true ]; then
    echo "--- Cleaning intermediate files ---" | tee -a "$LOG"
    rm -f "$TMP1_R1" "$TMP1_R2" "$TMP2_R1" "$TMP2_R2" \
          "$TMP3_R1" "$TMP3_R2" "$TMP4_R1" "$TMP4_R2"
fi

echo "[$(date)] Pipeline finished successfully." | tee -a "$LOG"
echo "Final Outputs:" | tee -a "$LOG"
echo "  R1: $R1_OUT" | tee -a "$LOG"
echo "  R2: $R2_OUT" | tee -a "$LOG"
```

保存并退出文本编辑器，运行脚本

```bash
bash 4_trim_data_bbduk.sh
```

对于sample CRR1071866而言，output结果（bbduk_CRR1071866_20251208_1327.log）如下：

| 步骤                               | 描述                                                | 实际操作                                       | 评价                                    |
| ---------------------------------- | --------------------------------------------------- | ---------------------------------------------- | --------------------------------------- |
| 1: ftm=5 (Force Trim Modulo 5)     | 强制修剪末端，使读长为 5 的倍数                     | 修剪了2.99% reads 的末端，移除 0.05% bases     | 轻微调整                                |
| 2: ftl=10 (Force Trim Left 10bp)   | 切除每个读的前 10bp（修正测序起始不稳定）           | 强制移除所有 reads 的左端 10bp，6.70% bases    |                                         |
| 3: Adapter Trimming (去除接头)     | 使用adapters参考去除Illumina等adapter               | 总移除 0.20% bases，82 reads                   | 序列低污染，有效处理 paired-end overlap |
| 4: PhiX Filtering (去除 PhiX 污染) | 使用 phix 参考过滤（k=31, hdist=1）                 | 0 reads/bases，即无污染                        | 0% 移除，库构建无 spike-in 残留         |
| 5: Quality Filtering (质量过滤)    | 修剪低质量端 (qtrim=rl, trimq=15)，丢弃 <30bp reads | Q-trim 6.00% reads，移除 0.11% bases，64 reads | 保留高质量数据                          |

#### 第五步：Bbduk_Trim_raw_data FastQC

对BBDuk 管道运行后序列进行FastQC质控

```bash
cd ../5_fastqc_quality_data
nano 5_fastqc_quality_data.sh
```

进入`5_fastqc_quality_data.sh`，输入以下script

```bash
#!/bin/bash
INPUT_DIR="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/2_qc/test/4_trim_data_bbduk"
fastqc -t 4 -f fastq -o . ${INPUT_DIR}/*.fastq
```

保存并退出文本编辑器，运行脚本

```bash
bash 5_fastqc_quality_data.sh
```

打开本地Terminal，下载并查看.html文件

```bash
scp u24149@logini.tongji.edu.cn:"/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/2_qc/test/5_fastqc_quality_data/*.html" /Users/xu.sj/Desktop
#输入服务器密码
```

![image-20251209103629231](/Users/xu.sj/Library/Application Support/typora-user-images/image-20251209103629231.png)

BBDuk处理后，所有碱基（A/C/G/T）从头到尾基本平直，无明显峰谷。

> [!NOTE]
>
> 由于空气样本主要捕获游离微生物DNA、病毒、细菌、古菌等，宿主细胞DNA很少进入空气滤膜，因此无需使用bowtie2做host removal。

### 正式运行

创建正式运行目录与样本列表

```bash
# 1. 切换到项目根目录
cd /ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1

# 2. 创建正式运行的目录结构 (与 test 平行)
mkdir -p 2_qc/formal
mkdir -p 3_megahit/formal
mkdir -p 4_mapping/formal
mkdir -p 5_binning/formal

# 3. 生成样本ID列表 (自动获取 CRR 开头的文件夹名)
ls 1_rawdata | grep "CRR" > sample_list.txt
```

将 FastQC、Fastp 和 BBDuk 步骤合并为一个脚本 `run_qc_batch.sh`

```bash
nano ./2_qc/run_qc_batch.sh
```

```bash
#!/bin/bash
#SBATCH -J QC_Batch
#SBATCH -p intel
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=10G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

SAMPLE_ID=$1
set -e

# ================= 路径配置 =================
BASE_DIR="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1"
RAW_DIR="${BASE_DIR}/1_rawdata/${SAMPLE_ID}"
OUT_DIR="${BASE_DIR}/2_qc/formal/${SAMPLE_ID}"

# 输入文件 (根据 Tree 结构)
R1_IN="${RAW_DIR}/${SAMPLE_ID}_f1.fq.gz"
R2_IN="${RAW_DIR}/${SAMPLE_ID}_r2.fq.gz"

mkdir -p ${OUT_DIR}/fastp ${OUT_DIR}/bbduk ${OUT_DIR}/fastqc_final

# ================= 1. Fastp (PolyG/Adapter) =================
source /share/home/u24149/miniconda3/etc/profile.d/conda.sh
conda activate fastp

echo "[${SAMPLE_ID}] Starting Fastp..."
fastp \
  -i "${R1_IN}" -I "${R2_IN}" \
  -o "${OUT_DIR}/fastp/${SAMPLE_ID}_fp_1.fq.gz" \
  -O "${OUT_DIR}/fastp/${SAMPLE_ID}_fp_2.fq.gz" \
  --detect_adapter_for_pe \
  --trim_poly_g --trim_poly_x \
  --cut_tail --cut_window_size 4 --cut_mean_quality 20 \
  --length_required 75 \
  --thread 8 \
  --html "${OUT_DIR}/fastp/${SAMPLE_ID}.fastp.html" \
  --json "${OUT_DIR}/fastp/${SAMPLE_ID}.fastp.json"

conda deactivate

# ================= 2. BBDuk (Deep Cleaning) =================
conda activate bbmap

# [还原] 保持您原本的变量设置 
ADAPTERS="adapters"
PHIX="phix"
MEM="8g"

echo "[${SAMPLE_ID}] Starting BBDuk..."

# 这里为了批量运行效率，将 Step 1-5 串联执行（逻辑不变），引用还原的变量
bbduk.sh -Xmx${MEM} \
  in="${OUT_DIR}/fastp/${SAMPLE_ID}_fp_1.fq.gz" \
  in2="${OUT_DIR}/fastp/${SAMPLE_ID}_fp_2.fq.gz" \
  out="${OUT_DIR}/bbduk/${SAMPLE_ID}_clean_1.fq" \
  out2="${OUT_DIR}/bbduk/${SAMPLE_ID}_clean_2.fq" \
  ref="${ADAPTERS},${PHIX}" \
  ftm=5 ftl=10 \
  ktrim=r k=23 mink=11 hdist=1 tbo tpe \
  qtrim=rl trimq=15 minlength=30 \
  overwrite=t t=8

# 压缩输出以节省空间
gzip -f "${OUT_DIR}/bbduk/${SAMPLE_ID}_clean_1.fq"
gzip -f "${OUT_DIR}/bbduk/${SAMPLE_ID}_clean_2.fq"

conda deactivate

# ================= 3. FastQC Final =================
echo "[${SAMPLE_ID}] Starting Final FastQC..."
fastqc -t 8 -o "${OUT_DIR}/fastqc_final" \
  "${OUT_DIR}/bbduk/${SAMPLE_ID}_clean_1.fq.gz" \
  "${OUT_DIR}/bbduk/${SAMPLE_ID}_clean_2.fq.gz"

echo "[${SAMPLE_ID}] QC Pipeline Finished."
```

```bash
#提交任务
cd 2_qc
mkdir -p logs
cat ../sample_list.txt | while read id; do sbatch run_qc_batch.sh $id; done
```



---

## 3.Megahit 

### 前期环境准备

- [ ] 已安装megahit——environment location：/share/home/u24149/miniconda3/envs/megahit

### Test: 以CRR1071866为例

```bash
cd /ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1 #切换到自己工作目录下
mkdir -p ./3_megahit/test 
cd ./3_megahit/test
nano run_assembly_megahit.sh
```

进入`run_assembly_megahit.sh`，输入以下script

```bash
#!/bin/bash
# MEGAHIT Assembly Script for CRR1071866 (Final Tuned Version)

# 遇到错误立即停止
set -e

# ================= 1. 环境激活 =================
source /share/home/u24149/miniconda3/etc/profile.d/conda.sh
conda activate megahit

# ================= 2. 输入路径 (Step 4 BBDuk 结果) =================
INPUT_DIR="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/2_qc/test/4_trim_data_bbduk"
R1="${INPUT_DIR}/CRR1071866_qc.R1.fastq"
R2="${INPUT_DIR}/CRR1071866_qc.R2.fastq"

# ================= 3. 输出路径 =================
TARGET_PARENT_DIR="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/3_megahit/test"

# 输出文件夹
OUTDIR="${TARGET_PARENT_DIR}/CRR1071866_megahit_out"

# ================= 4. 参数设置 =================
THREADS=16 # 线程数（匹配 SLURM CPU）
MEMORY=150000000000 # 150 GB（绝对值，单位字节）
MIN_LEN=500 # 过滤掉 <500bp 的序列

# ================= 5. 开始运行 =================
LOG="megahit_CRR1071866_$(date +%Y%m%d_%H%M).log"
echo "[$(date)] Starting MEGAHIT assembly..." | tee "$LOG"
echo "Output Dir: $OUTDIR" | tee -a "$LOG"

# 监控可用 RAM
echo "Available RAM before run: $(free -h --si | grep 'Mem:')" | tee -a "$LOG"

# 移除 shell 内存限制
ulimit -v unlimited

# 清理旧临时文件（防干扰）
rm -rf "$OUTDIR/tmp"

# 运行 MEGAHIT（用 --k-list 强制短列表，省内存）
megahit \
  -1 "$R1" \
  -2 "$R2" \
  -t "$THREADS" \
  -m "$MEMORY" \
  --min-contig-len "$MIN_LEN" \
  --k-list "21,41,61,81" \ #在基因组组装中，软件把你的测序数据打断成一个个小片段（K-mer）来拼图。小的 K值容易拼错，搞不定重复序列。大的 K值能跨过很长的重复序列，拼出来的 Contig 通常更长。但极其吃内存，而且如果测序深度不够，反而拼不出来。
  --out-prefix "CRR1071866" \
  -o "$OUTDIR" \
  -f 2>&1 | tee -a "$LOG"
 
 # ================= 6. 结果确认 =================
FINAL_FILE="${OUTDIR}/CRR1071866.contigs.fa"
if [ -f "$FINAL_FILE" ]; then
    echo "------------------------------------------------" | tee -a "$LOG"
    echo "[$(date)] Assembly Finished Successfully!" | tee -a "$LOG"
    echo "Final Result: $FINAL_FILE" | tee -a "$LOG"
   
    # 统计条数
    COUNT=$(grep -c ">" "$FINAL_FILE")
    echo "Total Contigs (>= 500bp): $COUNT" | tee -a "$LOG"
    
    # 额外统计（可选，用 seqkit 如果可用）
    if command -v seqkit &> /dev/null; then
        seqkit stats "$FINAL_FILE" | tee -a "$LOG"
    fi
else
    echo "[$(date)] Error: Output file not found! Check internal log: $OUTDIR/CRR1071866.log" | tee -a "$LOG"
fi
```

保存并退出，赋予执行权限

```bash
chmod +x run_assembly_megahit.sh
```

创建并提交 SLURM 作业（由于最开始直接bash失败，内存不足（Exit Code -9），所以切换到计算节点）

```bash
nano submit_megahit.sh
```

```bash
#!/bin/bash
#SBATCH --job-name=megahit_CRR1071866
#SBATCH --partition=fata  # 高内存首选；可改 amd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G  # 200 GB，覆盖峰值 + 裕度
#SBATCH --time=04:00:00
#SBATCH --output=slurm_megahit_%j.out
#SBATCH --error=slurm_megahit_%j.err

# 运行 MEGAHIT 脚本
bash run_assembly_megahit.sh
echo "Job $SLURM_JOB_ID finished at $(date)"
```

保存并退出，赋予执行权限并提交

```bash
chmod +x submit_megahit.sh
sbatch submit_megahit.sh
```

运行结束后可以在test/位置下tail -n 50 slurm_megahit_812622.out查看结束部分，以CRR1071866为例，Total Contigs (>= 500bp): 334960；由于畜牧场宏基因组多样性高，易产生多 contigs。该数值合理。

- #### MEGAHIT 核心结果如下： （位置：CRR1071866_megahit_out/）

| 文件                  | 含义                                                         | 用途                                    |
| --------------------- | ------------------------------------------------------------ | --------------------------------------- |
| CRR1071866.contigs.fa | 最终 contigs FASTA 文件（所有 k-mer 迭代合并后的连续序列）；内置过滤 ≥500 bp | 基因注释（Prokka）、丰度分析（Kraken2） |
| checkpoints.txt       | 组装检查点记录                                               | 追踪哪个阶段耗时长                      |
| CRR1071866.log        | MEGAHIT 内部详细日志（e.g., k-mer 统计、内存使用、错误）     | 检查 N50（组装质量指标）、内存峰值等等  |
| intermediate_contigs/ | 中间 contigs 子目录，按 k-mer 值（k21, k41, k61, k81）分文件夹。每个 k 值生成临时组装结果，用于迭代扩展 | 可删除（rm -rf），省空间                |
| done                  | 空标记文件，表示组装完全成功                                 | 用于脚本验证                            |

### 正式运行



---

## 4.Mapping

### 前期环境准备

- [ ] 已安装samtools——environment location：/share/home/u24149/miniconda3/envs/samtools
- [ ] 已安装metabat2——environment location：/share/home/u24149/miniconda3/envs/metabat2
- [ ] bbmap在第二步QC已安装

### Test: 以CRR1071866为例

#### 总体思路：

> 第一步将比对产生的乱序 BAM 文件进行排序（Sort）和建立索引（Index）;
>
> 第二步读取排序好的 BAM 文件，计算每条 Contig（基因组片段）上覆盖了多少 Reads（即 **丰度/深度**），生成一个 depth_bbmap.txt 文件。
>
> #MetaBAT2 分箱软件依靠“深度信息”来判断哪些 Contig 属于同一个细菌（原理是：同一个细菌的不同基因片段，在同一样本里的丰度接近）。

#### 第一步：Mapping, Sort, index

```bash
cd /ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/4_mapping #切换到自己工作目录下
mkdir ./test 
cd ./test
nano run_assembly_megahit.sh
```

```bash
#!/bin/bash
#SBATCH -J Map_CRR      # 任务名叫 Map_CRR
#SBATCH -p amd          # 指定投递到 amd 队列
#SBATCH -N 1            # 占用 1 个节点（比对软件不能跨节点，申请多了也是浪费）
#SBATCH -n 8            # 申请 8 个 CPU 核心
#SBATCH --mem=50G       # 申请 50GB 内存
#SBATCH -o %j.log       # 日志输出到ID.log
#SBATCH -e %j.err       # 错误日志输出到ID.err

# ==========================================
# Mapping Pipeline for CRR1071866
# 步骤：BBMap 比对 -> Samtools 排序 -> Samtools 索引
# 输出目录：.../4_mapping/test
# ==========================================

# 遇到错误立即停止
set -e
echo "[$(date)] Job Starting on node: $(hostname)"

# ================= 1. 定义路径  =================
OUTDIR="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/4_mapping/test"
# --- [输入1: 参考基因组] (来自 Step 3 Megahit) ---
REF="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/3_megahit/test/CRR1071866_megahit_out/CRR1071866.contigs.fa"
# --- [输入2: 测序 Reads] (来自 QC 第4步 BBDuk) ---
R1="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/2_qc/test/4_trim_data_bbduk/CRR1071866_qc.R1.fastq"
R2="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/2_qc/test/4_trim_data_bbduk/CRR1071866_qc.R2.fastq"
# --- [中间文件名定义] ---
RAW_BAM="${OUTDIR}/CRR1071866.bbmap.bam"
SORTED_BAM="${OUTDIR}/CRR1071866.sorted.bam"
STATS="${OUTDIR}/CRR1071866.stats.txt"

# 资源配置
# Java 内存设置要略小于申请的总内存 (50G)，留点给系统。如果占满了 50G，任务会被立刻 Kill 掉。
JAVA_MEM="45g"
THREADS=8

# ================= 2. 运行 BBMap (比对) =================
echo "------------------------------------------------"
echo "[$(date)] Step 1: Running BBMap..."

# 激活 bbmap 环境
source /share/home/u24149/miniconda3/etc/profile.d/conda.sh
conda activate bbmap

# 检查输入是否存在
if [ ! -f "$REF" ]; then
    echo "Error: Reference file not found: $REF"
    exit 1
fi

# 执行比对
bbmap.sh -Xmx$JAVA_MEM \ # 告诉 Java 虚拟机：最大可以用 45G 内存
  ref="$REF" \ # 参考基因组
  in="$R1" \ # 双端测序数据的左端
  in2="$R2" \ # 双端测序数据的右端
  out="$RAW_BAM" \ # 输出结果：RAW_BAM (未排序的 BAM 文件)
  minid=0.90 \ # 【重要】相似度阈值。只有匹配度 > 90% 才算比对成功，过滤掉乱比对的垃圾。
  statsfile="$STATS" \ # 生成一个统计报告（比如有多少 reads 比对上了）
  threads=$THREADS \
  overwrite=t # 如果文件已存在，直接覆盖，不要报错

echo "[$(date)] Mapping finished."

# ================= 3. 运行 Samtools (排序与索引) =================
echo "------------------------------------------------"
echo "[$(date)] Step 2: Running Samtools Sort & Index..."

# 切换环境
conda activate samtools

# 排序
echo "Sorting BAM..."
samtools sort -@ $THREADS -o "$SORTED_BAM" "$RAW_BAM"

# 索引
echo "Indexing BAM..."
samtools index "$SORTED_BAM"

# ================= 4. 清理与完成 =================
# 删除未排序的大文件以节省空间
rm -f "$RAW_BAM"

echo "------------------------------------------------"
echo "[$(date)] Pipeline Finished Successfully!"
echo "Final Result: $SORTED_BAM"
ls -lh "$SORTED_BAM" "$SORTED_BAM.bai"
```

保存并退出，提交任务

```bash
sbatch run_assembly_megahit.sh
```

任务结束后，首先看看日志文件最后有没有报错

```bash
tail -n 20 815206.log #更改为自己的作业号
# 成功标志：Pipeline Finished Successfully!
# 失败标志：Killed、Segmentation fault...
```

然后检查比对统计，`CRR1071866.stats.txt`会显示多少Reads成功比对上了

```bash
cat CRR1071866.stats.txt
```

结果&注释如下：

```bash
Reads Used:           	98276990	(13652618917 bases) #读长总数（总碱基数）

Mapping:          	2606.817 seconds. #映射时间
Reads/sec:       	37700.00
kBases/sec:      	5237.28

Pairing data:   	pct pairs	num pairs 	pct bases	   num bases

mated pairs:     	 58.2309% 	 28613783 	 58.2673% 	  7955011256
bad pairs:       	  4.5863% 	  2253625 	  4.6135% 	   629867934
insert size avg: 	  254.36 #平均插入大小为 254.36 bp，这在 Illumina 短读长库中是典型的
#评估：配对率 58% 属于中等水平（理想 >70%，取决于库构建和物种复杂度）。坏配对低，说明重复区域或污染不多。

Read 1 data:      	pct reads	num reads 	pct bases	   num bases  #正向读长统计

mapped:          	 69.8628% 	 34329520 	 69.9205% 	  4776070202 #映射率
unambiguous:     	 69.6298% 	 34215059 	 69.6892% 	  4760271844 #无歧义映射，高，表明特异性好
ambiguous:       	  0.2329% 	   114461 	  0.2313% 	    15798358 #歧义映射，极低，好
low-Q discards:  	  0.0000% 	        0 	  0.0000% 	           0
perfect best site:	 35.4250% 	 17407305 	 35.4384% 	  2420698400 
semiperfect site:	 36.4574% 	 17914612 	 36.4731% 	  2491371814
#完美/半完美位点：35.43% + 36.46% = 71.89%（大部分读长精确对齐）
rescued:         	  0.7223% 	   354942

Match Rate:      	      NA 	       NA 	 98.1500% 	  4706854670 #匹配率：98.15%，极高，序列保真度好
Error Rate:      	 47.9149% 	 16448971 	  1.5879% 	    76150927
Sub Rate:        	 47.6667% 	 16363745 	  1.1572% 	    55496663 #替换
Del Rate:        	  1.1225% 	   385344 	  0.4067% 	    19503659 #缺失
Ins Rate:        	  1.0040% 	   344682 	  0.0240% 	     1150605 #插入
N Rate:          	  2.5830% 	   886739 	  0.2621% 	    12568264
#评估：Read 1 映射优秀，错误主要为替换（常见于测序噪声），总体保真度高。

Read 2 data:      	pct reads	num reads 	pct bases	   num bases

mapped:          	 68.7809% 	 33797897 	 68.8368% 	  4695980159
unambiguous:     	 68.5593% 	 33689032 	 68.6168% 	  4680976026
ambiguous:       	  0.2215% 	   108865 	  0.2199% 	    15004133
low-Q discards:  	  0.0000% 	        0 	  0.0000% 	           0

perfect best site:	 24.2795% 	 11930573 	 24.2631% 	  1655206534
semiperfect site:	 25.0261% 	 12297437 	 25.0111% 	  1706233928
rescued:         	  0.8512% 	   418247

Match Rate:      	      NA 	       NA 	 97.7386% 	  4607665120
Error Rate:      	 63.6903% 	 21525968 	  2.0155% 	    95017454
Sub Rate:        	 63.5212% 	 21468831 	  1.6046% 	    75646955
Del Rate:        	  1.1008% 	   372058 	  0.3880% 	    18291298
Ins Rate:        	  0.9878% 	   333870 	  0.0229% 	     1079201
N Rate:          	  2.4846% 	   839743 	  0.2458% 	    11588883
#评估：Read 2 略逊于 Read 1，但错误率仍低，适合下游分析。
```

最后检查.bam文件大小，防止空文件

```bash
ls -lh CRR1071866.sorted.bam #结果显示文件有5.0G
#进一步检查核心output文件：CRR1071866.sorted.bam
source /share/home/u24149/miniconda3/etc/profile.d/conda.sh
conda activate samtools # 先激活 samtools 环境
samtools view -H CRR1071866.sorted.bam | head #前十行内容
conda deactivate # 退出环境
```

结果&注释如下：

```bash
@HD	VN:1.4	SO:coordinate 
#@HD: Header 头文件；VN:1.4: 格式版本号；SO:coordinate：Sort Order = Coordinate，明确提示这个文件已经被排过序了。只有带有这个标记，后续的 Binning 软件才肯认它。

@SQ	SN:k81_75160 flag=1 multi=3.0000 len=547	LN:547 #@SQ: Sequence
#SN:k81_75160: Sequence Name。MEGAHIT 组装出来的 Contig 名字（k81 代表它是用 k-mer 81 组装出来的）。通常在命令行里不加任何关于 K 值的参数时，MEGAHIT 内部的代码逻辑是这样的：对于标准的双端测序（通常是 150bp），它会自动设置 K-mer 列表为： 21, 29, 39, 59, 79, 99, 119, 141。由于我这里最开始没有用集群算，设置了81的参数，所以这里是k81。
#LN:547: Length。这条 Contig 长 547bp。

@SQ	SN:k81_450954 flag=1 multi=2.0000 len=672	LN:672
@SQ	SN:k81_977068 flag=1 multi=3.0000 len=733	LN:733
@SQ	SN:k81_1052226 flag=0 multi=1.1833 len=643	LN:643
@SQ	SN:k81_300636 flag=1 multi=1.0000 len=755	LN:755
@SQ	SN:k81_225477 flag=0 multi=7.0051 len=869	LN:869
@SQ	SN:k81_150319 flag=1 multi=2.0000 len=687	LN:687
@SQ	SN:k81_901909 flag=1 multi=1.0000 len=744	LN:744
@SQ	SN:k81_676431 flag=1 multi=3.0000 len=1101	LN:1101
```

#### 第二步：Depth_bbmap

```bash
cd /ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/4_mapping/test/
nano summarize_bam_depths.sh
```

```bash
#!/bin/bash
#SBATCH -J Depth_CRR1071866
#SBATCH -p amd
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=30G
#SBATCH -o %j.log
#SBATCH -e %j.err

# ==========================================
# Depth Calculation Pipeline for CRR1071866
# 步骤：读取 Sorted BAM -> 计算 Contig 丰度 -> 生成 depth.txt
# ==========================================

# 遇到错误立即停止
set -e
echo "[$(date)] Job Starting on node: $(hostname)"

# ================= 1. 定义路径 =================
# 工作目录：继续放在 4_mapping/test 下，或者你可以指向 5_binning 目录
WORK_DIR="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/4_mapping/test"

# --- [输入: 已排序的 BAM 文件] (来自 Mapping 步骤) ---
BAM_FILE="${WORK_DIR}/CRR1071866.sorted.bam"
ls
# --- [输出: 深度统计文件] ---
OUT_DEPTH="${WORK_DIR}/CRR1071866.depth.txt"

# ================= 2. 环境激活 =================
echo "------------------------------------------------"
echo "[$(date)] Activating Environment..."

source /share/home/u24149/miniconda3/etc/profile.d/conda.sh
conda activate metabat2

# ================= 3. 运行计算 (JGI Tool) =================
echo "------------------------------------------------"
echo "[$(date)] Step 1: Summarizing BAM contig depths..."

# 检查输入是否存在
if [ ! -f "$BAM_FILE" ]; then
    echo "Error: BAM file not found: $BAM_FILE"
    exit 1
fi

# --outputDepth: 指定输出文件名
jgi_summarize_bam_contig_depths \
    --outputDepth "$OUT_DEPTH" \
    "$BAM_FILE"

echo "[$(date)] Calculation finished."

# ================= 4. 结果确认 =================
echo "------------------------------------------------"
if [ -f "$OUT_DEPTH" ]; then
    echo "[$(date)] Pipeline Finished Successfully!"
    echo "Depth File: $OUT_DEPTH"
    # 打印前几行看看长什么样
    echo "Head of output file:"
    head -n 5 "$OUT_DEPTH"
else
    echo "Error: Output file was not generated."
    exit 1
fi
```

```bash
sbatch summarize_bam_depths.sh
```

查看output`CRR1071866.depth.txt`文件 

```bash
head -n 5 CRR1071866.depth.txt
```

```bash
contigName	contigLen	totalAvgDepth	CRR1071866.sorted	CRR1071866.sorted-var
k81_75160	547	3.00875	3.00875	1.74114
k81_450954	672	2.63402	2.63402	2.55946
k81_977068	733	3.99689	3.99689	5.52959
k81_1052226	643	7.32849	7.32849	24.13

#contigName: 基因组片段名称
#totalAvgDepth : 覆盖深度。比如 k81_1052226 的深度是 7.3，说明这个片段在测序中被覆盖了 7 次左右。
```

### 正式运行

----

## 5.Binning & Taxonomy

### 前期环境准备

- [ ] 已安装checkm2——environment location：/share/home/u24149/miniconda3/envs/checkm2
- [ ] 已安装gtdbtk——environment location：/share/home/u24149/miniconda3/envs/gtdbtk (需手动指定python版本)
- [ ] metabat2在第四步Mapping已安装
- [ ] CheckM2参考数据库手动下载，路径：/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/reference_database/checkm2_database/checkm2_v3/CheckM2_database/uniref100.KO.1.dmnd
- [ ] GTDB-Tk 参考数据库（25年4月最新版本，release226近140G）手动下载，路径：/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/reference_database/GTDB/release226

### Test: 以CRR1071866为例

#### 总体思路

> 分箱 (Binning) $\rightarrow$ CheckM2 质检脚本 $\rightarrow$ 物种注释 (Taxonomy)
>
> 第一步：根据 Mapping 的深度信息和序列特征，拆分成一个个独立的细菌Bins
>
> 第二步：对分出来的 Bins 进行质量打分，查看Completeness (完整度) 和 Contamination (污染度)
>
> 第三步：将高质量的 Bins 确定物种分类

#### 第一步：Binning

```bash
cd /ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/5_binning/test/
nano submit_metabat2_CRR1071866.sh
```

```bash
#!/bin/bash
#SBATCH -J CheckM2_CRR1071866
#SBATCH -p intel
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=20G
#SBATCH -o %j.log
#SBATCH -e %j.err

set -e

# ================= 1. 环境与路径 =================
source /share/home/u24149/miniconda3/etc/profile.d/conda.sh
conda activate metabat2

# 输入文件
CONTIGS="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/3_megahit/test/CRR1071866_megahit_out/CRR1071866.contigs.fa"
DEPTH="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/4_mapping/test/CRR1071866.depth.txt"

# 输出目录
OUT_DIR="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/5_binning/test/CRR1071866_metabat"
mkdir -p "$OUT_DIR"

# 输出前缀
PREFIX="${OUT_DIR}/CRR1071866.bin"

# ================= 2. 运行 MetaBAT2 =================
echo "[$(date)] Running MetaBAT2..."

# -m 1500: 保留大于1500bp的contig (默认是2500，设小一点能多捡回点数据)
metabat2 \
  -i "$CONTIGS" \
  -a "$DEPTH" \
  -o "$PREFIX" \
  -t 32 \
  -v

echo "[$(date)] MetaBAT2 finished!"
echo "Generated Bins:"
ls -lh "${PREFIX}"*.fa | head
```

```bash
sbatch submit_metabat2_CRR1071866.sh
```

查看结果：

```bash
cd /ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/5_binning/test/CRR1071866_metabat
ls -lh CRR1071866.bin*.fa | wc -l   # bins 数量
ls -lh CRR1071866.bin*.fa          # 列出每个 bin 的大小
```

#### 第二步：CheckM2质检

```bash
cd /ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/5_binning/test/
nano submit_checkm2_CRR1071866.sh
```

```bash
#!/bin/bash
#SBATCH -J Bin_CRR1071866
#SBATCH -p intel
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=20G
#SBATCH -o %j.log
#SBATCH -e %j.err

set -e
echo "[$(date)] Job Starting..."

# ================= 1. 环境与数据库 =================
source /share/home/u24149/miniconda3/etc/profile.d/conda.sh
conda activate checkm2

# 设置数据库路径
export CHECKM2DB="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/reference_database/checkm2_database/checkm2_v3/CheckM2_database/uniref100.KO.1.dmnd"

# ================= 2. 输入输出路径 =================
# 输入：指向 MetaBAT2 的结果文件夹
INPUT_DIR="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/5_binning/test/CRR1071866_metabat"
# 输出：CheckM2 结果文件夹
OUT_DIR="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/5_binning/test/CRR1071866_checkm2"
mkdir -p "$OUT_DIR"

# ================= 3. 运行 CheckM2 =================
echo "------------------------------------------------"
echo "[$(date)] Running CheckM2 Prediction..."
echo "Database: $CHECKM2DB"
echo "Input Bins: $INPUT_DIR"

# 检查输入文件夹是否为空
if [ -z "$(ls -A $INPUT_DIR/*.fa 2>/dev/null)" ]; then
    echo "Error: No .fa files found in $INPUT_DIR. Did MetaBAT2 finish successfully?"
    exit 1
fi

# 运行预测
checkm2 predict \
  --threads 16 \
  --input "${INPUT_DIR}"/*.fa \
  --output-directory "$OUT_DIR" \
  --force

echo "------------------------------------------------"
echo "[$(date)] CheckM2 Finished!"
echo "Report generated at: $OUT_DIR/quality_report.tsv"
```

```bash
sbatch submit_checkm2_CRR1071866.sh
```

查看结果：

```bash
cd /ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/5_binning/test/CRR1071866_checkm2

# 完整报告（用 column -t 对齐）
cat quality_report.tsv | column -t
```

判断标准：

| 质量                  | Completeness  | Contamination | 目标数量               | 评价                     |
| :-------------------- | :------------ | :------------ | :--------------------- | :----------------------- |
| High-quality          | ≥90%          | ≤5%           | ≥5-10 个最好，越多越好 | 优（可直接下游分析）     |
| Medium-quality        | ≥70%          | ≤10%          | 20-50 个常见           | 良（可用于物种水平分析） |
| Medium-quality (宽松) | ≥50%          | ≤10%          | -                      | 可接受                   |
| Low-quality           | <50% 或高污染 | -             | 大量                   | 差（需优化 binning）     |

快速筛选命令：

```bash
# High-quality MAGs
awk 'NR==1 || ($2 >= 90 && $3 <= 5)' quality_report.tsv | column -t

# Medium + High（>=70% Comp, <=10% Cont，常用阈值）
awk 'NR==1 || ($2 >= 70 && $3 <= 10)' quality_report.tsv | column -t
```

#### 第三步：GTDB-Tk 物种注释

```bash
cd /ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/5_binning/test/
nano submit_gtdbtk_CRR1071866.sh
```

```bash
#!/bin/bash
#SBATCH -J GTDB_CRR
#SBATCH -p intel
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mem=250G
#SBATCH -o %j.log
#SBATCH -e %j.err

set -euo pipefail
echo "[$(date)] Job Starting..."

# ================= 1. 环境与数据库 =================
source /share/home/u24149/miniconda3/etc/profile.d/conda.sh
conda activate gtdbtk

# GTDB 数据库的文件夹路径
export GTDBTK_DATA_PATH="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/reference_database/GTDB/release226"

# ================= 2. 输入输出路径 =================
BASE_DIR="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1"

# 输入：依然是 MetaBAT2 的结果 (或者你可以写脚本只挑出 CheckM2 质量好的 bins)
INPUT_DIR="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/5_binning/test/CRR1071866_metabat"

# 输出
OUT_DIR="/ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/5_binning/test/CRR1071866_gtdbtk"
mkdir -p "$OUT_DIR"

# ================= 3. 运行 GTDB-Tk =================
echo "------------------------------------------------"
echo "[$(date)] Checking GTDB Database..."
if [ -z "$GTDBTK_DATA_PATH" ]; then
    echo "Error: You explicitly explicitly set the GTDBTK_DATA_PATH variable!"
    exit 1
fi
echo "Database Path: $GTDBTK_DATA_PATH"

echo "[$(date)] Running GTDB-Tk Classify Workflow..."

# classify_wf 是全自动流程 (Ani -> Align -> Classify)
gtdbtk classify_wf \
    --genome_dir "$INPUT_DIR" \
    --extension fa \
    --out_dir "$OUT_DIR" \
    --cpus 16 \
    --skip_ani_screen \
    --pplacer_cpus 1  # 限制 pplacer 线程以防内存爆炸

echo "------------------------------------------------"
echo "[$(date)] GTDB-Tk Finished!"
echo "Results saved in: $OUT_DIR"
```

```bash
sbatch submit_gtdbtk_CRR1071866.sh
```

查看结果：

```bash
cd /ssdfs/datahome/u24149/Lab_Member/Students/Shengjie.Xu/Metagenomics/Livestock_farm_1/5_binning/test/CRR1071866_gtdbtk

# 统计有多少个 Bin 失败了
wc -l identify/gtdbtk.failed_genomes.tsv

# 统计有多少个 Bin 成功注释到了“种”水平 (含有 s__ 且后面有名字)
grep "s__" gtdbtk.bac120.summary.tsv | grep -v "s__;" | wc -l

# 提取分类结果中的“界”和“门”，然后统计数量
cut -f 2 gtdbtk.bac120.summary.tsv | cut -d ';' -f 1,2 | sort | uniq -c | sort -nr
```

### 正式运行

----

> [!TIP]
>
> 友情分享：Typora 下载link🔗:https://dn25g9a610.feishu.cn/drive/folder/Qjkif8yasljLind2L5GcKNDAnjg
