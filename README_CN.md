# AKmerGWAS

一套基于 k-mer 的全基因组关联分析（GWAS）流水线，可接受任意输入格式——短读长序列、组装结果或完整基因组。

该流水线使用 [kmtricks](https://github.com/tlemane/kmtricks) 统计 k-mer，
将存在/缺失矩阵转换为 PLINK bed 分块，对每个分块和性状运行
[GEMMA](https://github.com/genetics-statistics/GEMMA) 关联分析，并用
`bowtie` 或 `bwa` 将显著的 k-mer 比对回参考基因组。

**1.0.0** 版是将原先的 Python/Snakemake/Bash 技术栈彻底重写为单一纯 Python 编排器（仅依赖标准库），**完全由命令行驱动**——不存在配置文件。参见
[从 0.1.0 迁移](#migration-from-010)。

---

## 安装

```bash
conda install -c bioconda akgwas
akgwas --version
```

运行时外部工具（以 conda 依赖形式声明）：
`kmtricks`、`gemma`、`plink`、`bowtie`、`bwa`。

---

## 快速开始

该流水线会将显著的 k-mer 比对回参考基因组，因此需要一个**参考索引**。akgwas 从不自行构建索引——在默认比对器（`bowtie`）下，需为每个参考基因组用 `bowtie-build` 构建一次，并将生成的索引**前缀**写入 `ref_list.txt`：

```bash
# 1. build a bowtie index (one-off; reuse across runs)
mkdir -p /data/indexes
bowtie-build /data/genomes/pseudo10001.fasta /data/indexes/pseudo10001
# -> writes /data/indexes/pseudo10001.{1,2,3,4}.ebwt and .rev.{1,2}.ebwt

# 2. ref_list.txt: one index prefix per line -- the path you passed to
#    bowtie-build, WITHOUT the .1.ebwt / .rev.1.ebwt suffix
echo /data/indexes/pseudo10001 > /data/run/ref_list.txt
```

（若使用 `--aligner bwa`，则改为列出 bwa 索引前缀——参见 [输入](#inputs) 下的 `ref_list.txt` 小节。）

然后，将 `fof.txt`、`pheno.txt` 和 `ref_list.txt` 放入 `--programhome`，并运行整条流水线：

```bash
# Required files in --programhome: fof.txt, pheno.txt, ref_list.txt
akgwas all --programhome /data/run --cores 32
```

这就是全部接口。每个值都有默认值；传入对应的命令行参数（flag）即可覆盖它。只有 `--programhome`（相对路径和所有输出的基准目录）是必需的。

单独运行各个阶段：

```bash
akgwas build     --programhome /data/run --cores 8    # kmtricks count/merge + bed chunking
akgwas kgwas     --programhome /data/run --cores 32   # K/PCA + GEMMA + collect
akgwas alignment --programhome /data/run --cores 8    # map significant k-mers
```

试运行（不实际执行任何操作），或在忽略已保存进度的情况下重新运行：

```bash
akgwas all --programhome /data/run --dry-run
akgwas all --programhome /data/run --force
```

对某个子命令运行 `--help` 可查看分组后的选项（输入/输出、build、multibed、qk、kgwas、alignment、collect），例如
`akgwas build --programhome /data/run --block-size 5G --maf 0.05`。

---

## 输入

三个纯文本文件位于 `--programhome` 中（文件名可用 `--fof`、
`--pheno`、`--ref-list` 覆盖）。

### 1. `fof.txt` — 序列

每行一个样本：`name : file1 ; file2 ; ...`。

```
A88  : /data/pseudo88_1.fasta.gz ; /data/pseudo88_2.fasta.gz
A108 : /data/pseudo108.fasta.gz
```

### 2. `pheno.txt` — 表型表

单张表。**第一列存放样本 ID**，**表头行存放性状名称**，缺失值写为 `NA`。

```
ID      ft10     ft16
A88     61.25    22.1
A108    70       NA
A139    68.75    19.7
```

对于 GWAS，该流水线会：验证样本 ID 的*集合*与 bed 的 `.fam` 集合一致，将行重排为 `.fam` 的顺序，并为每个性状导出一个单列、无表头的文件（`<run_root>/pheno_single/<trait>.txt`）供 GEMMA `-p` 使用。

### 3. `ref_list.txt` — 参考索引

每行一个条目。在 `--aligner bowtie`（默认）下，条目为 **bowtie 索引前缀**——即你传给 `bowtie-build <ref.fasta> <prefix>` 的 `<prefix>`。
akgwas **不会**运行 `bowtie-build`；请先准备好索引，再将其路径写到这里：

```bash
bowtie-build /data/genomes/pseudo10001.fasta /data/indexes/pseudo10001
```

```
# ref_list.txt
/data/indexes/pseudo10001
/data/indexes/pseudo10002
```

在 `--aligner bwa` 下，条目为 bwa 索引前缀（或 `*.sa` 路径）。

---

## 命令行即配置

每个配置值都是一个带默认值的命令行参数（flag）；不从任何文件中读取配置。确切的调用方式会记录在运行清单（manifest）中，因此运行保持可复现。要点如下：

| 命令行参数（flag） | 默认值 | 含义 |
|---|---|---|
| `--programhome` | *(必需)* | 路径和输出的基准目录 |
| `--outdir` | `programhome` | 输出根目录（相对路径 → programhome） |
| `--run-name` | *(无)* | 按表型隔离输出，见下文 |
| `--kmer-len` | `31` | k-mer 长度 |
| `--hard-min` | `2` | 判定为存在的最小样本内计数（组装结果/基因组使用 `1`） |
| `--block-size` / `--min-block-size` / `--max-block-size` | `1G` / `16M` / `10G` | 分块大小边界，建议不小于500M |
| `--maf` | `0.01` | 次等位基因频率（MAF）过滤 |
| `--bed-backend` | `python` | `python`（直接生成 bed）或 `plink` |
| `--use-external-qk` `--kinship` `--pca` | 关闭 | 提供外部 K/PCA |
| `--sample-fraction` `--seed` `--gk` | `0.2` / `123456789` / `1` | 内部 K/PCA |
| `--threshold` | `0.0001` | 用于提取的 `p_wald` 阈值 |
| `--aligner` | `bowtie` | `bowtie` 或 `bwa` |
| `--max-hits` | `5000` | 丢弃命中数超过该值的 k-mer |
| `--gemma-combine` | 关闭 | 同时写出全部 k-mer 的完整转储（`gemma_combine/`；可能非常大） |

运行 `akgwas <subcommand> --help` 可查看完整列表。

### 块大小 —— `1G` / `16M` / `10G`

`--block-size` 及其上下界 `--min-block-size` / `--max-block-size` 是 `multibed`
阶段的**字节预算**。该阶段把 kmtricks 的文本矩阵流式切成一块块 PLINK 文件
（`bed/matrix_XXXX.{bed,bim,fam}`），每一块就是一个 GEMMA 的处理单位。所以块
大小是在“峰值磁盘/内存”和“GEMMA 调用次数”之间取舍：块越小 → 块数越多 → 峰值
占用越低，但跑的次数越多。

后缀是 **1024 进制（二进制）**，与 GNU `parallel --block` 一致
（`K = 1024`、`M = 1024²`、`G = 1024³`、`T = 1024⁴`），因此默认值对应：

| 参数 | 默认 | 字节数 | 约合 |
|---|---|---|---|
| `--block-size` | `1G` | 1,073,741,824 | 1 GiB ≈ 1.07 GB |
| `--min-block-size` | `16M` | 16,777,216 | 16 MiB ≈ 16.8 MB |
| `--max-block-size` | `10G` | 10,737,418,240 | 10 GiB ≈ 10.74 GB |

`--block-size` 是**目标值**：累计字节达到该预算时，当前块封口、开下一个块。取值
必须落在 `[--min-block-size, --max-block-size]` 之间。所有尺寸参数都用同一套后缀
（例如 `--block-size 5G`）。

> 注：当前实际切块只由 `--block-size` 决定；`--min-block-size` /
> `--max-block-size` 仅用于校验取值范围（`min ≤ block ≤ max`），未参与实际切块。

### 跨表型运行复用 bed/QK — `--run-name`

`build`、`multibed` 和 `qk` 是**与性状无关**的，因此它们的输出（`bed/`、`QK/`、`kmatrixs_out/`）直接位于 `outdir` 下并被共享。所有依赖表型的内容（`pheno_single/`、`temp_gemma/`、
`kmergwas_outputs/`，以及 kgwas/collect/alignment 的完成标记）位于 `run_root` = `outdir`（无运行名）或 `outdir/<run-name>` 下。

这样一来，你就可以把开销高昂的矩阵 + K/PCA 只计算一次，然后在不同时间运行多个表型集合，而不会相互覆盖，也不会错误地跳过阶段：

```bash
# once: matrix + K/PCA + first trait set
akgwas all  --programhome /data/run --cores 32

# later: a different phenotype table, reusing bed/ and QK/
akgwas kgwas     --programhome /data/run --run-name ft16 --pheno pheno_ft16.txt
akgwas alignment --programhome /data/run --run-name ft16
```

`ft16` 的输出落在 `/data/run/ft16/kmergwas_outputs/` 下，拥有自己的完成标记文件，因此不会与第一次运行冲突（也不会因第一次运行而被跳过）。

---

## 样本对齐——唯一至关重要的不变式

一切都以 **bed 的 `.fam` 顺序**为准，而该顺序本身源自
`fof.txt`（样本声明顺序——*并非*排序后的顺序）：

* 表型表被重排为 `.fam` 顺序（基于 ID 匹配）；
* 亲缘关系（kinship）/PCA 矩阵必须按照**与最终 GWAS bed 分块相同的样本顺序**构建；
* 工具会强制执行这一点：内部 k-mer K/PCA 通过拼接 bed 分块构建（绝不使用 `plink --merge-list`，因为它会按字典序重排样本），并且 PCA 的 `eigenvec` ID 会与规范顺序进行核对。

**当你提供外部 K/PCA 时**（`--use-external-qk`），你需要自行负责样本顺序。如果你的 PCA 文件带有样本 ID，请加上 `--pca-has-id`，工具会对其进行重排/验证；亲缘关系矩阵（没有 ID）则假定遵循相同的顺序——参见 `docs/troubleshooting.md`。

<span style="color:#e5484d">**最佳的保证此问题不会出错的方式是：将fof的样本名称按照字典顺序排序，之后表型和kinship、PCA的样本顺序也按照字典顺序，这样可以保证样本不会错配**</span>


---

## 输出

```
<outdir>/
├── bed/                     matrix_0001.{bed,bim,fam} ... (trait-independent)
├── QK/                      QK_kinship.txt, QK_pca1-3.txt, QK_samples.txt
├── .akgwas/                 akgwas.log, run_manifest[_<run>].json, <stage>.done
└── <run_root>/              <run_root> = . (no --run-name) or <run-name>/
    ├── pheno_single/        <trait>.txt (single column, fed to GEMMA)
    └── kmergwas_outputs/
        ├── sigk/            <trait>_sigK.txt          (k-mer af beta se p_wald)
        ├── siggeno/         <trait>_bed.{bed,bim,fam} (PLINK)
        ├── fasta/           <trait>_fasta.txt         (sequences for alignment)
        ├── gemma_combine/   <trait>.assoc_combined.txt  (only with --gemma-combine)
        ├── logs/            <trait>.log_combined.txt
        ├── samfile/         <trait>_<ref>.sam
        └── sigk_mapped/     <trait>_<ref>_chrpos_sigk.txt
```

`build` 阶段只保留 kmtricks 运行目录和 bed 分块：数 GB 的文本矩阵以及所有 `.tped` 中间文件从不写入磁盘（以流式 / 一次一个分块的方式处理）。

---

## 资源、磁盘与可复现性

* **磁盘。** 主要开销是 kmtricks 运行目录加上最终的 bed。
  分块大小（`--block-size`）限定了瞬时的峰值磁盘占用。
* **内存 / CPU。** k-mer 统计受分块大小约束；GEMMA 每个作业单线程运行，并发度来自 `--cores`。
* **可复现性。** 分块抽样种子为 `--seed`；工具版本、完整解析后的配置以及所选工具的版本会被写入
  `.akgwas/run_manifest[_<run-name>].json`。
* **续跑。** 每个阶段都会写入一个原子性的 `<stage>.done` 完成标记（带有输入指纹）。重新运行时，输入未发生变化的阶段会被跳过；`--force` 则忽略这些完成标记。依赖性状的阶段按 `--run-name` 分别标记。

---

## 从 0.1.0 迁移

| 方面 | 0.1.0 | 1.0.0 |
|---|---|---|
| 配置 | `config.yaml` | **无——使用命令行参数（flag）** |
| 表型 | `pheno_dir/` 下每个性状一个单列文件 | 一个多列 `pheno.txt`（ID + 性状列） |
| 编排 | Snakemake（`.smk`）+ GNU parallel/rush + Bash | 纯 Python |
| 默认比对器 | `bwa mem` | `bowtie -a -v 0`（也可用 bwa `-k31 -T31 -a -c 5000`） |
| CLI | `--build/--kgwas/...` 参数 | `akgwas build\|kgwas\|alignment\|all` |
| 共享 vs 每次运行 | 所有输出在单一目录树中 | `bed/`+`QK/` 共享；按性状的输出位于 `--run-name` 下 |
| 移除的依赖 | snakemake、perl、rush、parallel | — |

---

## 开发 / 测试

```bash
pip install -e .
pytest tests/              # unit tests run without external tools
pytest tests/test_smoke.py # end-to-end, needs kmtricks/gemma/plink/bowtie
python -m tests.make_synthetic /tmp/syn   # inspect a generated dataset
```
