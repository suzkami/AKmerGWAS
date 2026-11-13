# AKmerGWAS 运行示例（v6）

本文给出可直接照抄的命令：一条 `akgwas all` 跑完整流程、分三步单独运行（共用同一份 bed）、以及使用外部 bed / 外部 K/PCA 时的写法，最后是完整参数清单和每步必需输入。

> 相关文档：`docs/AKmerGWAS_Q&A.md`（外部 bed 的硬条件、v6 为什么不能直接跳过 PCA 等）。
> 文中的默认值与行为描述对应 v6 源码；不同版本可能变化。

---

## 1. 准备输入

示例目录结构：

```text
/data/akgwas_demo/
├── fof.txt                 # 样本 -> 输入文件列表
├── pheno.txt               # 表型表（多列）
├── ref_list.txt            # 参考基因组索引列表
├── reads/                  # 原始 reads / 组装 / 基因组，任意格式
│   ├── A88.fq.gz
│   └── A108.fq.gz
└── refs/                   # bowtie 或 bwa 索引前缀
    ├── TAIR10.1.bt2
    ├── TAIR10.2.bt2
    └── ...
```

`fof.txt`：每行一个样本，格式 `样本名 : 文件1 ; 文件2 ; ...`，`#` 开头为注释。

```text
A88  : /data/akgwas_demo/reads/A88.fq.gz
A108 : /data/akgwas_demo/reads/A108.fq.gz
# 双端示例：S1 : /path/S1_R1.fq.gz ; /path/S1_R2.fq.gz
```

`pheno.txt`：第一列样本 ID，第一行是表头（性状名），缺失值写 `NA`。

```text
ID      ft10    ft16
A88     61.25   22.1
A108    70      NA
```

`ref_list.txt`：每行一个索引前缀。bowtie 用 `bowtie-build` 的前缀；bwa 用 `bwa index` 的前缀（bwa 模式下条目以 `.sa` 结尾也可，程序会自动去掉）。

```text
/data/akgwas_demo/refs/TAIR10
```

必须保证三处样本标识一致：`fof.txt` 的样本名、`pheno.txt` 的 ID 列、bed `.fam` 的 IID 列。`--programhome` 必须是绝对路径，其他相对路径都相对它解析。

---

## 2. 一条命令跑完整流程：`akgwas all`

`all` 依次执行 build → multibed → qk → kgwas → collect → alignment。

```bash
akgwas all \
    --programhome /data/akgwas_demo \
    --cores 32 \
    --kmer-len 31 \
    --fof fof.txt \
    --pheno pheno.txt \
    --hard-min 2 \
    --recurrence-min 5 \
    --partitions 0 \
    --build-threads 8 \
    --multibed-threads 8 \
    --block-size 1G \
    --maf 0.01 \
    --bed-backend python \
    --sample-fraction 0.2 \
    --seed 123456789 \
    --gk 1 \
    --qk-threads 4 \
    --threshold 1e-4 \
    --ref-list ref_list.txt \
    --aligner bowtie \
    --alignment-threads 4 \
    --gemma-combine
```

要点：

- `--cores` 只控制两个并行池：kgwas 的 (chunk × trait) 任务、alignment 的参考基因组任务。build / multibed / qk 各自的并行度由 `--build-threads`、`--multibed-threads`、`--qk-threads` 控制。
- 想先看看会跑什么而不真正执行：加 `--dry-run`。
- 中途失败后重跑同一条命令，已完成阶段会按 `.akgwas/` 里的完成标记跳过；要强制重跑加 `--force`。
- `--gemma-combine` 会把全部 k-mer 的关联结果也导出（`gemma_combine/`，可能非常大），默认关闭。

---

## 3. 分三步运行（共用同一份 bed）

三步必须使用相同的 `--programhome`（以及相同 `--outdir` / `--bed-dir` / `--run-name`，如果改过）。build 专有参数在 kgwas / alignment 中传入不会报错，但不会生效。

### 第 1 步：`akgwas build`（count + 切块）

```bash
akgwas build \
    --programhome /data/akgwas_demo \
    --cores 32 \
    --kmer-len 31 \
    --fof fof.txt \
    --hard-min 2 \
    --recurrence-min 5 \
    --partitions 0 \
    --build-threads 8 \
    --multibed-threads 8 \
    --block-size 1G \
    --maf 0.01 \
    --bed-backend python
```

产出 `<programhome>/kmatrixs_out/` 和 `<programhome>/bed/matrix_0001.{bed,bim,fam}` …，以及 `bed/sample_id_inorder.txt`。这里的 `matrix_*` 就是后面 kgwas 与 alignment 共用的 bed。

### 第 2 步：`akgwas kgwas`（QK + 关联 + 汇总）

```bash
akgwas kgwas \
    --programhome /data/akgwas_demo \
    --cores 32 \
    --fof fof.txt \
    --pheno pheno.txt \
    --sample-fraction 0.2 \
    --seed 123456789 \
    --gk 1 \
    --qk-threads 4 \
    --threshold 1e-4 \
    --gemma-combine
```

这一步包含三个内部阶段：qk（用上一步的 bed 抽样建 `QK/QK_kinship.txt` 与 `QK/QK_pca1-3.txt`）、kgwas（每个 chunk × 每个性状跑一次 GEMMA）、collect（汇总 sigK、fasta、siggeno 等）。

### 第 3 步：`akgwas alignment`（回贴参考基因组）

```bash
akgwas alignment \
    --programhome /data/akgwas_demo \
    --cores 8 \
    --fof fof.txt \
    --kmer-len 31 \
    --ref-list ref_list.txt \
    --aligner bowtie \
    --alignment-threads 4
```

必须注意：

- `--kmer-len` 必须与第 1 步相同（默认都是 31）。不一致会导致 CIGAR `<k>M` / `NM:i:0` 过滤把全部命中丢弃，输出为空。
- v6 的 preflight 会无条件检查 fof 文件存在，所以 alignment 也要带上 `--fof`；内容不参与 alignment，但文件必须在。
- 多次 GWAS 可共用 bed 与 QK：用 `--run-name ft10`、`--run-name ft16` 分别跑 kgwas，输出隔离在 `<outdir>/<run-name>/`；alignment 要用同一个 `--run-name`。

---

## 4. 外部 bed / 外部 K/PCA 的写法

### 4.1 只用外部 K/PCA（bed 仍由 akgwas 生成）

`akgwas all` 可以直接使用外部 K/PCA：

```bash
akgwas all \
    --programhome /data/akgwas_demo \
    --cores 32 \
    --fof fof.txt \
    --pheno pheno.txt \
    --recurrence-min 5 \
    --maf 0.01 \
    --use-external-qk \
    --kinship K.txt \
    --pca PCA.txt \
    --pca-has-id \
    --threshold 1e-4 \
    --ref-list ref_list.txt \
    --aligner bowtie
```

分步运行时，只在 kgwas 步加外部参数：

```bash
akgwas build --programhome /data/akgwas_demo --cores 32 \
    --fof fof.txt --recurrence-min 5 --maf 0.01

akgwas kgwas --programhome /data/akgwas_demo --cores 32 \
    --fof fof.txt --pheno pheno.txt \
    --use-external-qk --kinship K.txt --pca PCA.txt --pca-has-id

akgwas alignment --programhome /data/akgwas_demo --cores 8 \
    --fof fof.txt --kmer-len 31 --ref-list ref_list.txt
```

文件格式：

- **K（亲缘矩阵）**：GEMMA `-k` 的格式，首行样本数 n，随后 n 行 × n 列数值。行/列顺序必须与 bed fam 顺序一致，v6 不校验 K 的样本 ID。
- **PCA / 协变量**：加 `--pca-has-id` 时按第 2 列做样本 ID 匹配，取第 3–5 列作为协变量（最多 3 列）。不加时视为无表头、每行若干数值列，只校验行数等于样本数，行序按 bed fam。

### 4.2 外部 bed + 外部 K/PCA（不能跑 `all`）

把外部 bed 放到 `<outdir>/bed/` 后，**不要运行 `build`、`multibed` 或 `all`**：multibed 会删除 `bed/matrix_*` 并重建。等价的“下游 all”是 kgwas 模块 + alignment：

```bash
akgwas kgwas \
    --programhome /data/akgwas_demo \
    --cores 32 \
    --fof fof.txt \
    --pheno pheno.txt \
    --use-external-qk \
    --kinship K.txt \
    --pca PCA.txt \
    --pca-has-id \
    --threshold 1e-4

akgwas alignment \
    --programhome /data/akgwas_demo \
    --cores 8 \
    --fof fof.txt \
    --kmer-len 31 \
    --ref-list ref_list.txt \
    --aligner bowtie
```

如果希望用一条命令串起来，可以这样写（`&&` 保证 kgwas 成功后再比对）：

```bash
akgwas kgwas --programhome /data/akgwas_demo --cores 32 \
    --fof fof.txt --pheno pheno.txt \
    --use-external-qk --kinship K.txt --pca PCA.txt --pca-has-id \
  && akgwas alignment --programhome /data/akgwas_demo --cores 8 \
    --fof fof.txt --kmer-len 31 --ref-list ref_list.txt
```

外部 bed 的硬条件（详见 `docs/AKmerGWAS_Q&A.md`）：

1. 文件为 `<outdir>/bed/matrix_XXXX.bed/.bim/.fam`，PLINK1 SNP-major（`plink --make-bed` 产物即可）。
2. 所有 chunk 的 `.fam` 样本顺序完全一致；建议写 `bed/sample_id_inorder.txt`（每行一个 IID，按 fam 顺序）。
3. **`.bim` 第 2 列必须是 k-mer 序列本身**，否则 collect 取不到变异，`siggeno` 与 alignment 结果为空。
4. `fof.txt` 必须存在，样本名与 `pheno.txt`、fam IID 一致（不跑 build 时，fof 里的文件路径可写占位）。
5. 如果沿用旧 outdir，加 `--force` 或换新 outdir，避免完成标记让 qk / kgwas 直接跳过。

### 4.3 外部 bed + 内部 QK

不传 `--use-external-qk`，qk 会从外部 bed 抽样重建 K/PCA：

```bash
akgwas kgwas \
    --programhome /data/akgwas_demo \
    --cores 32 \
    --fof fof.txt \
    --pheno pheno.txt \
    --sample-fraction 0.2 \
    --seed 123456789 \
    --gk 1 \
    --qk-threads 4

akgwas alignment \
    --programhome /data/akgwas_demo \
    --cores 8 \
    --fof fof.txt \
    --kmer-len 31 \
    --ref-list ref_list.txt
```

---

## 5. 全部参数清单

说明：四个子命令（`build` / `kgwas` / `alignment` / `all`）共用同一套参数定义。**参数层面唯一必填的是 `--programhome`**，其余不传就使用默认值。“生效步骤”表示该参数在哪一步真正起作用；传给其他步骤不会报错，但会被忽略。

### 5.1 通用 / 输入输出

| 参数 | 默认值 | 生效步骤 | 必填/可选 | 说明 |
| --- | --- | --- | --- | --- |
| `--programhome` | 无 | 全部 | **必填** | 绝对路径，所有相对路径与输出的基准目录 |
| `--cores` | `1` | kgwas、alignment | 可选 | 并行任务数：kgwas 的 (chunk × trait)，alignment 的参考基因组数 |
| `--dry-run` | 关闭 | 全部 | 可选 | 只打印将要执行的命令，不写任何文件 |
| `--force` | 关闭 | 全部 | 可选 | 忽略完成标记，强制重跑 |
| `-q`, `--quiet` | 关闭 | 全部 | 可选 | 减少屏幕输出（错误仍会打印） |
| `--outdir` | 等于 programhome | 全部 | 可选 | 输出根目录；相对路径相对 programhome |
| `--run-name` | 空 | kgwas、collect、alignment | 可选 | 把表型相关输出隔离到 `<outdir>/<run-name>/`，bed 与 QK 仍共享 |
| `--fof` | `fof.txt` | build、kgwas；alignment 仅检查存在 | 可选 | 样本列表文件 |
| `--pheno` | `pheno.txt` | kgwas | 可选 | 多列表型文件 |
| `--kmer-len` | `31` | build、alignment | 可选 | k-mer 长度；alignment 必须与建 bed 时一致 |
| `--ulimit-min` | `4096` | 全部 preflight | 可选 | 要求的最小 open-file limit |
| `--bed-dir` | `bed` | build、qk、kgwas、collect | 可选 | 共享的 bed chunk 目录名 |
| `--qk-dir` | `QK` | qk、kgwas | 可选 | 共享的 K/PCA 目录名 |
| `--temp-dir` | `temp_gemma` | kgwas、collect | 可选 | 每次运行下 GEMMA 中间结果目录名 |
| `--final-dir` | `kmergwas_outputs` | collect、alignment | 可选 | 每次运行的最终结果目录名 |
| `--kmtricks-run-dir` | `kmatrixs_out` | build、multibed | 可选 | kmtricks run 目录名 |

### 5.2 build（k-mer 计数）

| 参数 | 默认值 | 生效步骤 | 必填/可选 | 说明 |
| --- | --- | --- | --- | --- |
| `--hard-min` | `2` | build | 可选 | 样本内最小计数，达到才算“存在” |
| `--recurrence-min` | `2` | build | 可选 | 跨样本最小重复数，同时传给 kmtricks 的 count 与 merge |
| `--partitions` | `0` | build | 可选 | kmtricks 分区数；`0` 表示自动 |
| `--build-threads` | `8` | build | 可选 | kmtricks count 线程数（merge 固定 2 线程） |
| `--aggregate-to-stdout`, `--no-aggregate-to-stdout` | 开启 | build / multibed | 可选 | multibed 直接消费 aggregate 的 stdout；关闭则先写临时文件再读 |

### 5.3 multibed（切块，属于 build 模块）

| 参数 | 默认值 | 生效步骤 | 必填/可选 | 说明 |
| --- | --- | --- | --- | --- |
| `--block-size` | `1G` | multibed | 可选 | 每个 chunk 的目标字节预算，单位按 1024 进制（K/M/G/T） |
| `--min-block-size` | `16M` | multibed | 可选 | 只参与参数校验（要求 min ≤ block-size ≤ max），运行时不会截断或合并 chunk |
| `--max-block-size` | `10G` | multibed | 可选 | 同上 |
| `--maf` | `0.01` | multibed | 可选 | 次等位基因频率过滤：`max(#absent, #present) <= int(n - maf*n)` 才保留 |
| `--multibed-threads` | `8` | multibed | 可选 | `kmtricks aggregate` 的线程数 |
| `--bed-backend` | `python` | multibed | 可选 | `python` 直接写 bed；`plink` 先写 tped 再由 plink 转换 |

### 5.4 qk（亲缘矩阵与 PCA，属于 kgwas 模块）

| 参数 | 默认值 | 生效步骤 | 必填/可选 | 说明 |
| --- | --- | --- | --- | --- |
| `--use-external-qk` | 关闭 | qk | 可选 | 使用外部 K/PCA；开启时必须同时提供 `--kinship` 与 `--pca` |
| `--kinship` | 空 | qk | 条件必填 | 外部亲缘矩阵路径（相对 programhome） |
| `--pca` | 空 | qk | 条件必填 | 外部 PCA / 协变量文件路径（相对 programhome） |
| `--sample-fraction` | `0.2` | qk | 可选 | 内部 K/PCA 抽样使用的 chunk 比例 |
| `--seed` | `123456789` | qk | 可选 | 抽样随机种子（保证可复现） |
| `--gk` | `1` | qk | 可选 | GEMMA `-gk`：`1` 居中，`2` 标准化 |
| `--pca-has-id` | 关闭 | qk | 可选 | 外部 PCA 带样本 ID（第 2 列），按 ID 校验并重排 |
| `--qk-threads` | `1` | qk | 可选 | `plink --pca` 线程数 |

### 5.5 kgwas（关联分析）

| 参数 | 默认值 | 生效步骤 | 必填/可选 | 说明 |
| --- | --- | --- | --- | --- |
| `--threshold` | `0.0001` | kgwas | 可选 | 显著阈值：`p_wald <= threshold` 的 k-mer 写入 sigK |

### 5.6 collect（汇总）

| 参数 | 默认值 | 生效步骤 | 必填/可选 | 说明 |
| --- | --- | --- | --- | --- |
| `--gemma-combine`, `--no-gemma-combine` | 关闭 | collect | 可选 | 额外写出全部 k-mer 的关联汇总（`gemma_combine/`，可能很大） |

### 5.7 alignment（回贴参考基因组）

| 参数 | 默认值 | 生效步骤 | 必填/可选 | 说明 |
| --- | --- | --- | --- | --- |
| `--ref-list` | `ref_list.txt` | alignment | 可选 | 参考索引列表（相对 programhome） |
| `--aligner` | `bowtie` | alignment | 可选 | `bowtie` 或 `bwa` |
| `--max-hits` | `5000` | alignment | 可选 | 单个 k-mer 命中数超过该值时整条丢弃 |
| `--alignment-threads` | `1` | alignment | 可选 | 每个比对进程的线程数（bowtie `-p` / bwa `-t`） |

### 5.8 顶层

| 参数 | 说明 |
| --- | --- |
| `akgwas --version` | 打印版本 |
| `akgwas --errors` | 打印错误码表 |
| `akgwas --help` / `akgwas <子命令> --help` | 帮助 |

---

## 6. 每步必需输入速查

“必需输入的参数”严格来说只有 `--programhome`；下表的“必需数据”指该步骤运行前必须存在的文件或工具，可以放在默认路径，也可以用对应参数指向别处。

| 步骤 | 必需输入的参数 | 必须存在的数据 / 工具 | 产出 |
| --- | --- | --- | --- |
| `build` | `--programhome` | `--fof` 指定的文件（里面的输入文件也要存在）；`kmtricks`；`--bed-backend plink` 时还要 `plink` | `kmatrixs_out/`、`bed/matrix_*`、`bed/sample_id_inorder.txt` |
| `kgwas` | `--programhome` | `--fof`、`--pheno`；`bed/matrix_*.bed`（来自 build 或外部）；`gemma` 必需，内部 QK 另需 `plink`；外部 QK 另需 `--kinship` 与 `--pca` 文件 | `QK/`、`<run_root>/temp_gemma/`、`<run_root>/kmergwas_outputs/` |
| `alignment` | `--programhome` | `--fof` 文件（仅 preflight 检查存在）；`--ref-list`；上一层产出的 `kmergwas_outputs/sigk/*_sigK.txt` 与 `fasta/`；`bowtie` 或 `bwa` | `kmergwas_outputs/samfile/`、`kmergwas_outputs/sigk_mapped/` |
| `all` | `--programhome` | 以上全部（`fof`、`pheno`、`ref_list`，以及 kmtricks、gemma、plink、比对器） | 以上全部 |

---

## 7. 重跑、跳过与 `--force`

- 每个阶段成功后会在 `<outdir>/.akgwas/` 写完成标记；再次运行同一命令会跳过已完成阶段。
- 完成标记按输入文件的时间戳/大小判断，**不包含参数签名**。修改了影响结果的参数（如 `--hard-min`、`--recurrence-min`、`--maf`、`--kmer-len`）后，必须加 `--force`，必要时先删除 `kmatrixs_out/`、`bed/`、`QK/` 再重跑，避免新旧结果混用。
- 只想重跑某一步：单独执行对应子命令并加 `--force`。

---

## 8. 输出位置速查

| 路径 | 内容 |
| --- | --- |
| `<outdir>/kmatrixs_out/` | kmtricks 计数结果 |
| `<outdir>/bed/matrix_XXXX.{bed,bim,fam}` | 供 kgwas / alignment 使用的分块 bed |
| `<outdir>/bed/sample_id_inorder.txt` | 样本顺序 |
| `<outdir>/QK/QK_kinship.txt` | 亲缘矩阵（稳定命名） |
| `<outdir>/QK/QK_pca1-3.txt` | 前 3 个 PCA（稳定命名） |
| `<outdir>/QK/QK_samples.txt` | qk 使用的样本顺序 |
| `<run_root>/pheno_single/` | 每个性状一个单列表型文件 |
| `<run_root>/temp_gemma/` | 每个 chunk 的 GEMMA 中间结果 |
| `<run_root>/kmergwas_outputs/sigk/` | 每个性状的显著 k-mer（`<trait>_sigK.txt`） |
| `<run_root>/kmergwas_outputs/fasta/` | 显著 k-mer 的 FASTA，供比对 |
| `<run_root>/kmergwas_outputs/siggeno/` | 显著 k-mer 的基因型 bed |
| `<run_root>/kmergwas_outputs/samfile/` | 比对 SAM |
| `<run_root>/kmergwas_outputs/sigk_mapped/` | 带参考位置的结果（`<trait>_<ref>_chrpos_sigk.txt`） |
| `<outdir>/.akgwas/` | 完成标记、日志、run manifest |

`<run_root>` 为 `<outdir>`，或在指定 `--run-name` 时为 `<outdir>/<run-name>`。

---

## 9. 常见问题

- **`E601` / 找不到 `QK_kinship.cXX.txt`**：GEMMA 的真实输出名里 `XX` 是字面量。v6 已按 `cXX` / `sXX` 识别；若仍失败，检查 `QK/output/` 下的实际文件名与 GEMMA 版本。
- **`no bed chunks under .../bed`**：`--outdir` / `--bed-dir` 与上一步不一致，或外部 bed 命名不是 `matrix_*.bed`。
- **样本 ID 不一致**：fof、pheno、bed fam 三处必须完全一致；外部 bed 场景优先检查 fam 的 IID 列。
- **alignment 输出为空**：`--kmer-len` 与建库时不一致，或显著 k-mer 本身为空（可先看 `sigk_mapped/` 是否有文件、`samfile/` 是否为空）。
- **外部 bed 跑完没有 siggeno**：`.bim` 第 2 列不是 k-mer 序列（详见 `docs/AKmerGWAS_Q&A.md`）。
