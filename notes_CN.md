# AKmerGWAS v6 — 使用与排障笔记（note）

本文件记录 v6 的操作要点：`build` 失败后的自救步骤、几条常见用法的澄清
（QK 在哪个阶段建 / `--run-name` 隔离 / bowtie 与 bwa 各自比对 / `--kmer-len`），
以及命令行用法上的一条限制（子命令不能叠加）。

---

## 1. build 不成时如何重跑

如果 `build` 没成、或产物状态可疑，**删除 `kmatrixs_out/`、`bed/`、`QK/` 三个目录，
然后加 `--force` 重跑**：

```bash
rm -rf /data/run/kmatrixs_out /data/run/bed /data/run/QK
akgwas build --programhome /data/run --cores 8 --force
```

说明：

- `build` 阶段在调 `kmtricks pipeline` 前本身会 `remove_if_exists(kmatrixs_out/)`
  （kmtricks 拒绝已存在的 `--run-dir`），所以 `kmatrixs_out/` 通常能自清；
  真正需要手动删的是**残留的 `bed/` 和 `QK/`**。
- `--force` 的作用是**忽略完成标记、强制重跑各阶段**。
- 为什么要“删产物 + `--force`”一起做：完成标记只对**输入**取指纹，**不看产物是否
  存在、也不看参数是否变化**。若 `fof.txt` 没变而产物已被删，阶段会误判
  `valid -> skip`，随后撞上不存在的产物报错（典型：`build: marker valid -> skip`
  紧跟 `[E401] kmtricks aggregate ... [--run-dir .../kmatrixs_out] Directory doesn't
  exist!`）。删掉产物再加 `--force`，能把这种状态干净重置。

---

## 2. 常见用法澄清

### 2.1 QK 在哪个阶段创建；单独跑 build 会不会建 QK？

阶段划分（`pipeline.py` 的 `MODULE_STAGES`）：

| 模块 | 阶段 |
|---|---|
| `build` | `build` + `multibed` |
| `kgwas` | `qk` + `kgwas` + `collect` |
| `alignment` | `alignment` |
| `all` | 全部六个 |

- `akgwas build` = kmtricks 计数 + merge（生成 `kmatrixs_out/`）**再** multibed
  （`kmtricks aggregate` 流式切块 → `bed/matrix_XXXX.{bed,bim,fam}`）。
  **它不建 QK。**
- **QK 在 `qk` 阶段创建**，而 `qk` 是 `kgwas`（也是 `all`）的**第一个**阶段。
  产物：`QK/QK_kinship.txt`、`QK/QK_pca1-3.txt`、`QK/QK_samples.txt`
  （中间体 `QK_bed.*` 用完即删）。
- QK 与性状无关、是**共享**产物（放 `outdir/QK`），其完成标记也不按 run-name 分，
  所以**只建一次**，后续换 `--run-name` 直接复用。

### 2.2 每个表型如何起一个新的 run？目录在哪？temp_gemma / kmergwas_outputs 怎么处理？

用 `--run-name` 隔离：

```bash
akgwas kgwas     --programhome /data/run --run-name ft16 --pheno pheno_ft16.txt --cores 32
akgwas alignment --programhome /data/run --run-name ft16 --cores 8
```

- `--pheno` 相对 `programhome` 解析（也可给绝对路径）。
- `run_root = outdir/<run-name>`（`outdir` 默认就等于 `programhome`），因此本次运行
  目录是 **`/data/run/ft16/`**，其下三样：
  - `pheno_single/` — 每个 trait 一个单列文件，喂 GEMMA
  - `temp_gemma/` — 每个 chunk 一个子目录：`output/<trait>.assoc.txt`、`sigk/`、`fasta/`
  - `kmergwas_outputs/` — collect + alignment 的最终产物
- 共享目录 `kmatrixs_out/`、`bed/`、`QK/` 留在 `/data/run/` 顶层，不重算。
- 跑第二个表型时：
  - `--run-name` **不同** → `temp_gemma`、`kmergwas_outputs` 在各自的 run_root 下，
    完全隔离（推荐用法）。
  - **复用同一个 run-name**、只换 `--pheno` → 共用同一 `temp_gemma` /
    `kmergwas_outputs`；kgwas 会重跑（pheno 路径进了输入指纹，指纹变了），文件按名
    覆盖（原子替换），但**旧表型遗留、新表型没有的 trait 文件不会被自动清理**。
    所以**每个表型建议单独给 `--run-name`**。

### 2.3 同一表型分别做 bowtie / bwa 比对；只改 ref_list.txt 吗？ksize 固定吗？

比对输出都在 `<run_root>/kmergwas_outputs/` 下：

```
samfile/<trait>_<refname>.sam
sigk_mapped/<trait>_<refname>_chrpos_sigk.txt
```

`<refname>` = `ref_list.txt` 里那条 index 路径的 basename，去掉序列后缀
（`.fasta/.fa/.fasta.gz/.gz/.sa` …）；basename 重复时会加 `_1`、`_2`。

- **文件名里没有 aligner 标识**。若 bowtie 索引和 bwa 索引的 basename 相同
  （例如都叫 `pseudo10001`），两次比对的 `samfile`/`sigk_mapped` 会**同名互相覆盖**。
  想同时保留两套结果：用**不同的 `--run-name`**（如 `ft16_bowtie` / `ft16_bwa`），
  或把两个索引前缀命名成不同的 basename。
- 换成 bwa **不只是改 `ref_list.txt`**：还要加 `--aligner bwa`，并把 `ref_list.txt`
  的内容换成 **bwa 的索引**（bwa index 前缀，或 `*.sa` 路径）。
- `alignment` 的完成标记输入里**不含 aligner**，只含 `ref_list` 文件指纹。若改了
  `--aligner` 但 `ref_list.txt` 内容没变，阶段可能被 `valid -> skip`，此时要补
  `--force`；正常换了索引路径、文件内容变了，指纹会变，会自动重跑。
- **ksize 不固定、也不从矩阵里自动探测**，它就是 CLI 的 `--kmer-len`（默认 31）。
  `build` 用它当 `kmtricks --kmer-size`；`alignment` 用它两处：精确匹配要的 CIGAR
  长度 `<k>M`，以及 bwa 的 `-k <k> -T <k>`。因此**如果 build 时不是默认值，单独跑
  alignment 时必须传同样的 `--kmer-len`**；否则（例如 build 用 25、alignment 用默认
  31）全部命中会被过滤、`sigk_mapped` 为空，bwa 的 `-k/-T` 也会错。
  （`kgwas`/`collect`/`qk`/`multibed` 都不用 `kmer_len`。）

---

## 3. 子命令不能叠加

一个命令只能有**一个**子命令（argparse subparsers）。

- `akgwas kgwas alignment ...` → 报 `unrecognized arguments: alignment`。
- `akgwas --programhome /data/run ...`（不带子命令）→ 只打印帮助，什么都不做。
- 想一次跑完 kgwas 再 alignment，三种办法：
  1. 分两条命令；
  2. shell 串联：
     ```bash
     akgwas kgwas --programhome /data/run --run-name ft16 --pheno pheno_ft16.txt --cores 32 \
       && akgwas alignment --programhome /data/run --run-name ft16 --cores 8
     ```
  3. 用 `all`（跑 `build→multibed→qk→kgwas→collect→alignment` 六个阶段；已完成的
     `build`/`multibed`/`qk` 会被完成标记跳过，所以实际只多花 kgwas+collect+alignment）：
     ```bash
     akgwas all --programhome /data/run --run-name ft16 --pheno pheno_ft16.txt --cores 32
     ```
- 想让 kgwas 和 alignment 用**不同** `--cores`，就分两条命令（`all` 只能统一给一个
  `--cores`）。
