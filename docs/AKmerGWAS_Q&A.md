# AKmerGWAS 问答记录（v6）

> 适用版本：v6（纯 Python 重写版）。文中行号以 `v6/akgwas/` 下的文件为准，改动后可能偏移。
> 本文只记录问答结论，未改动任何代码。

---

## Q1：v6 里不想使用 PCA 或任何协变量，怎么办？直接空着 `--pca` 可以吗？

**结论：不行。** v6 目前没有“关闭 PCA / 协变量”的开关，`--pca` 留空不会关掉协变量，只会导致报错。能跳过的只有协变量 `-c`；亲缘矩阵 `-k` 是 `-lmm` 必需的，不能跳过。

### 为什么空着 `--pca` 会失败

| 位置 | 行为 |
| --- | --- |
| `config.py:46-52, 225-231` | `--pca` 只写入 `qk.pca`；`validate()`（`config.py:269-293`）没有“允许为空”的分支 |
| `qk.py:222-230` | `--use-external-qk` 会同时要求 `kinship` 与 `pca` 文件存在；`pca` 为空时 `os.path.exists("")` 为 False，报 `E_INPUT_MISSING: external PCA not found:` |
| `kgwas.py:34-46` | 无论内部还是外部，都固定读取 `QK/QK_pca1-3.txt`，缺一个就报 `E_INPUT_MISSING` |
| `kgwas.py:102-104` | GEMMA 命令行里的 `-c <pca>` 是硬编码的，没有条件分支 |

另外，内部 qk 的 PC 数也是硬编码的：`qk.py:123-127` 固定 `plink --pca 3`，所以“只保留 1–2 个 PC”同样需要改代码。

### 四种可选做法

1. **继续使用 PCA（推荐，成本最低）。**
   内部模式只抽 `--sample-fraction`（默认 20%）的 chunk 拼成 `QK_bed`，再跑 `plink --pca 3` 与 `gemma -gk`，对 k-mer 矩阵来说开销很小。外部模式用 `--use-external-qk --kinship K --pca PCA` 直接复用你信得过的矩阵。

2. **换成你自己的协变量（不是“不要”，而是“用别的”）。**

   ```bash
   akgwas kgwas --programhome /abs/run \
       --use-external-qk --kinship K.txt --pca covariates.txt --pca-has-id
   ```

   `covariates.txt` 就是 GEMMA `-c` 的格式：每行一个样本、若干数值列；v6 只是把它叫 PCA。带 `--pca-has-id` 时按第 2 列做样本 ID 匹配，取第 3–5 列作为协变量（最多 3 列，`qk.py:252-273`）；不带 ID 时只校验行数，顺序默认与 bed fam 一致（`qk.py:236-247`）。

3. **真正不用协变量：需要改代码（改动很小）。**
   - `kgwas.py`：把 pca 从“必须存在”改为可选，只在文件存在时给 GEMMA 加 `-c`。
   - `qk.py`：`--use-external-qk` 时允许只给 `--kinship`；内部模式增加跳过 PCA 的分支（或新增 `--no-pca`）。
   - `cli.py` / `config.py`：新增开关并更新 `validate()`，避免空路径被当成缺失文件。
   这属于行为变更，需要确认后再动手。

4. **不改代码的 hack：全 0 的假 PCA（不推荐用于正式结果）。**
   造一个 N 行 × 3 列、所有值为 0 的文件作为 `--pca` 传入，可以绕过存在性检查。但 GEMMA 仍会把它当协变量矩阵，全 0 列与截距共线，不同版本表现不一（可能报奇异/共线，或给出异常数值）。如果只想先跑通，可先用小数据确认 GEMMA 不报错，并与真正无协变量的结果对比后再用。

### 关于“也不用亲缘矩阵”

`-lmm 1` 必须有 `-k`，v6 不支持跳过 K。若想换 K，用 `--use-external-qk --kinship <你的矩阵>`。注意 K 的行/列顺序默认按 bed fam 顺序，v6 不校验 K 的样本 ID（PCA 有 `--pca-has-id`，K 没有对应开关）。

### 一点建议

即使 K 已控制亲缘关系，k-mer 数据里仍常有更粗的群体分层，完全不给协变量容易得到大量假阳性。若确实想去掉 PCA，建议先保留 1–2 个 PC 做对照，确认显著性结果不是由分层驱动的。

---

## Q2：使用外部构建好、已切块的 bed，能否在此基础上继续跑 kgwas 和 alignment？

**结论：可以。** `kgwas` 模块（qk → kgwas → collect）和 `alignment` 都不读 `kmatrixs_out/`，只从 `<outdir>/bed/` 按 `matrix_*.bed` 找 chunk。但有 7 条硬条件，其中第 3、4 条不满足时不会报错，而是“跑成功但结果为空”。

### 数据流：哪些步骤真正读 bed

- **qk**：在 `bed/` 找 chunk，抽样拼接成 `QK_bed`，跑 `plink --pca` 与 `gemma -gk`（或直接使用外部 K/PCA）。
- **kgwas**：glob `bed/matrix_*.bed`；用第一个 chunk 的 `.fam` 校验并重排表型；每个 (chunk, trait) 跑一次 GEMMA。
- **collect**：用 sigK 里的 k-mer 字符串去对应 chunk 的 `.bim` 取子集，生成 `siggeno`。
- **alignment**：不碰 bed，只使用 `kmergwas_outputs/sigk/*_sigK.txt`、`fasta/` 与 `ref_list.txt`。

### 硬性条件

1. **位置与命名。** `matrix_XXXX.bed/.bim/.fam` 放在 `<outdir>/bed/`（`outdir` 默认等于 `--programhome`，可用 `--outdir` 修改；目录名可用 `--bed-dir` 修改）。必须是 PLINK1 SNP-major bed（`plink --make-bed` 的产物即可）。fam 行建议为 `0 <IID> 0 0 0 0`（与 `multibed.py:75` 一致）。
2. **所有 chunk 的 `.fam` 样本顺序完全一致。** qk 会把选中的 chunk 直接按行拼接，顺序不一致报 `E_SAMPLE_ORDER_MISMATCH`（`qk.py:114-120`）。建议写 `bed/sample_id_inorder.txt`（每行一个 IID，按 fam 顺序），qk 会优先使用它（`qk.py:70-78`）；不写则退回使用第一个 chunk 的 fam。
3. **`.bim` 第 2 列（variant ID）必须是 k-mer 序列本身。** 这是最关键的一条。collect 建 `siggeno` 时用 k-mer 字符串匹配 bim ID（`bed.py:191`），GEMMA 的 `rs` 列也来自这个 ID。若 ID 是 `rs123`、`kmer_0001` 之类的编号，GEMMA 照跑、`sigK` 也有内容，但 `siggeno` 一条都取不到，alignment 也就没有可映射项——全程不报错。bim 需保持标准 6 列（`kgwas.py:24` 的列下标 `[1,6,7,8,11]` 依赖 GEMMA 标准输出）。
4. **fof 仍然必须提供。** preflight 无条件检查 `--fof` 文件存在（`pipeline.py:99`），且 kgwas/collect 会核对 fof 样本名与表型样本 ID 的数量和集合（`pipeline.py:127-142`）。不跑 build 时文件路径可以是占位，但格式必须是 `name : somefile`。
5. **表型文件的样本 ID 必须对上第一个 chunk 的 fam。** kgwas 用第一个 chunk 的 fam 做重排与导出（`kgwas.py:53-54`，`pheno.py`），对不上会直接报错。
6. **K/PCA 二选一。**
   - 推荐 `--use-external-qk --kinship K --pca PCA`（加 `--pca-has-id` 会按 ID 校验/重排）。
   - 不给外部时，qk 会从你的 chunk 按 `--sample-fraction` / `--seed` 抽样重建 `QK/QK_kinship.txt` 和 `QK/QK_pca1-3.txt`；`kgwas` 只认这两个稳定文件名。
7. **alignment 依赖 collect 的产物。** 需要 `kmergwas_outputs/sigk/*_sigK.txt` 与 `fasta/`，以及 `ref_list.txt`；`--kmer-len` 必须等于生成这批 k-mer 的 k（否则 CIGAR `<k>M` 与 `NM:i:0` 过滤会把命中全部丢掉）。若用过 `--run-name`，alignment 必须带同一个 `--run-name`。

### 运行命令

```bash
akgwas kgwas --programhome /abs/run --outdir /abs/run \
    --fof fof.txt --pheno pheno.txt \
    --use-external-qk --kinship K.txt --pca PCA.txt --cores 8

akgwas alignment --programhome /abs/run --outdir /abs/run \
    --ref-list ref_list.txt --kmer-len 31 --cores 8
```

### 两个必须避开的坑

- **不要跑 `build`、`multibed` 或 `all`。** `multibed` 开头会删除 `bed/matrix_*` 并全部重建（`multibed.py:70-71`）；`build` 会删除 `kmatrixs_out` 并要求 kmtricks。
- **注意完成标记会跳过。** 同一 outdir 之前跑过的话，`.akgwas/` 里的 marker 可能让 qk/kgwas 直接输出 “marker valid -> skip”。换新 outdir，或加 `--force`。

### 运行前自检（只读）

```bash
ls bed/matrix_*.bed | head
head -2 bed/matrix_0001.bim                 # 第 2 列应是 ACGT...，长度 = --kmer-len
for f in bed/*.fam; do cut -d' ' -f2 "$f" | md5sum; done | sort -u | wc -l   # 应为 1
```

### 附：bim 里不是 k-mer 怎么办

如果有 ID→k-mer 对照表，可对每个 chunk 批量改名后再跑，例如 `plink --bfile <chunk> --update-name <map> --make-bed --out <new>`（列顺序以 plink 文档为准），或直接重写 `.bim` 第 2 列（其余列保持不变）。改名前确认对照表覆盖所有变异，否则 collect 仍会漏掉未命中的 k-mer。

---

## Q3：`--cores` 与 `--alignment-threads` 是什么关系？能同时做到 kgwas 32 并发、alignment 8 个参考吗？

**结论：** `--cores` 不是全局线程数，它只决定两个任务池的宽度（kgwas 的 chunk × trait 任务数、alignment 的参考并行数）。两个阶段共用同一个 `--cores`，因此“kgwas 32 并发 + alignment 只跑 8 个参考”**无法用一条 `akgwas all` 做到**，必须把 kgwas 和 alignment 分开跑，各自给自己的 `--cores`。

### 各阶段的并发到底由谁控制

| 阶段 | 并发/线程由谁控制 | `--cores` 是否生效 | 默认值 |
| --- | --- | --- | --- |
| build（count） | `--build-threads` → `kmtricks -t` | 不生效 | 8 |
| build（merge） | 硬编码 `-t 2` | 不生效 | 2 |
| multibed（aggregate） | `--multibed-threads` → `kmtricks aggregate -t` | 不生效 | 8 |
| multibed（MAF/切块/写 bed） | 单线程 Python 循环 | 不生效 | 1 |
| qk（PCA） | `--qk-threads` → `plink --pca --threads` | 不生效 | 1 |
| qk（kinship） | `gemma -gk`，无线程参数 | 不生效 | 1 |
| kgwas | `ThreadPoolExecutor(max_workers=cores)`，每个 GEMMA 设 `OMP_NUM_THREADS=1` 等 | **生效** | 1 |
| collect | 单线程顺序遍历 trait / chunk | 不生效 | 1 |
| alignment | 池宽 = `min(cores, 参考数)`；每个比对进程再用 `--alignment-threads` | **生效（仅池宽）** | 1 |

### kgwas：`--cores` 近似等于总线程数

kgwas 的任务数是 chunk 数 × 性状数，每个任务是一个**单线程** GEMMA（环境里强制 `OMP_NUM_THREADS=1`、`OPENBLAS_NUM_THREADS=1`、`MKL_NUM_THREADS=1`），kgwas 没有 per-GEMMA 的线程开关。所以同时跑的 GEMMA 数 ≈ `min(--cores, chunk 数 × 性状数)`，每个占 1 个 CPU。

### alignment：`--cores` × `--alignment-threads` = 峰值线程数（重点）

- 池宽：`workers = min(--cores, 参考数)`，一个 worker 负责一个参考。
- 参考内部**按 trait 串行**（trait 循环是顺序的），所以同一个参考同时最多一个比对进程。
- 每个比对进程的线程数由 `--alignment-threads` 决定：bowtie 传给 `-p`，bwa 传给 `-t`。
- 因此峰值线程数 ≈ `min(--cores, 参考数) × --alignment-threads`，**代码里没有任何上限或超订阅保护**。

举例：`--cores 32`、参考数 ≥ 32、`--alignment-threads 4`，最多会同时存在 **32 个比对进程 × 4 线程 = 128 个可运行线程**。在 32 核机器上就是 4 倍超订阅，整体会变慢。

内存不要按 128 倍估算，而要按“同时存在的进程数”估算：

- 每个比对进程持有：一份参考索引 + 很小的一份显著 k-mer FASTA + SAM 输出缓冲；
- **bowtie** 的索引是 mmap 的，同一份索引文件在多个进程间共享物理页；但 32 个**不同**参考就是 32 份索引同时驻留；
- **bwa mem** 每个进程各自加载索引，32 个进程就是 32 份索引内存。参考大时（人类基因组 bwa 索引 5–8 GB 量级）确实可能 OOM；
- 线程数本身不线性放大内存（线程栈基本是虚拟内存），真正占内存的是进程数和索引大小。

### 想要 kgwas 32、alignment 8 个参考：分两步

```bash
# kgwas 模块：qk + kgwas + collect，32 个并发单线程 GEMMA
akgwas kgwas --programhome /data/akgwas_demo --cores 32 \
    --fof fof.txt --pheno pheno.txt

# alignment：池宽 8，即最多 8 个参考并行
akgwas alignment --programhome /data/akgwas_demo --cores 8 \
    --fof fof.txt --kmer-len 31 --ref-list ref_list.txt \
    --alignment-threads 4
```

此时 alignment 峰值线程 ≈ `min(8, 参考数) × 4 = 32`，内存则是最多 8 份参考索引同时在跑。

如果写成 `akgwas all --cores 32`，alignment 会拿到 `min(32, 参考数)` 的池宽，不会自动变成 8。

分步时注意：

- 两步使用相同的 `--programhome`（以及相同的 `--outdir` / `--bed-dir` / `--run-name`，如果改过）。
- alignment 仍然要求 `--fof` 文件存在（preflight 会检查，内容不参与比对）。
- 如果先跑过 `akgwas all`，再想用新的 `--cores` 重跑 alignment，必须加 `--force`，否则完成标记会让 alignment 直接 “marker valid -> skip”。
- v6 没有“按阶段设置 cores”的参数，需要不同并发度时就别用 `all`，用分步运行。

### 32 核机器上的参考配置

| 阶段 | 建议 |
| --- | --- |
| kgwas | `--cores 32`（每个 GEMMA 单线程） |
| build（count） | `--build-threads 32`（或 16，看磁盘） |
| build（merge） | 固定 2 线程，无需设置 |
| multibed（aggregate） | `--multibed-threads 8–16`；切块循环本身单线程 |
| qk | `--qk-threads 1–4`（PCA/kinship 本身并行度有限） |
| alignment（参考 ≥ 32） | `--cores 32 --alignment-threads 1`，或 `--cores 8 --alignment-threads 4` |
| alignment（只有 1 个参考） | `--cores` 无效（池宽 = 1），只设 `--alignment-threads 8–16` |
| alignment（参考大、内存紧张） | 降低 `--cores`（2–4），或分批跑 `ref_list.txt` |

一句话记忆：`--cores` 是“同时跑几个任务”，`--*-threads` 是“每个任务几线程”；kgwas 峰值线程 ≈ `min(cores, chunk × trait)`，alignment 峰值线程 ≈ `min(cores, 参考数) × alignment-threads`。调参时让峰值线程接近物理核数即可，多参考场景还要单独按“同时在跑的索引份数”核对内存。
