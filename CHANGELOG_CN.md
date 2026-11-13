# 更新日志

此处记录 AKmerGWAS 的所有重要变更。格式遵循
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/)，项目使用
[语义化版本](https://semver.org/)。

## [1.0.0] — 2026-10-12

将 0.1.0 的 Python/Snakemake/Bash 技术栈彻底重写为单一纯 Python 编排器，
**无任何第三方运行时依赖**（仅使用标准库）。外部工具（kmtricks、gemma、plink、bowtie、bwa）
仍以子进程方式调用。

### 新增

- 四个子命令：`akgwas build`、`akgwas kgwas`、`akgwas alignment`、
  `akgwas all`，以及 `--dry-run`、`--force`、`--version`、`--errors`。
- 每个配置项都对应一个命令行参数（flag），且均有默认值；只有
  `--programhome` 是必填项。`akgwas <subcommand> --help` 按功能分区列出全部参数。
- `--run-name`，用于把依赖表型的输出隔离到 `outdir/<name>/` 下，
  这样一份 `bed/`+`QK/`（只需计算一次）即可服务多个表型的运行，而不会
  相互覆盖，也不会错误地跳过阶段。
- 集中、稳定的错误码（`E1xx`–`E9xx`），配有易懂的提示，并有
  `docs/troubleshooting.md` 作为参考。
- 原子化的完成标记（`<stage>.done` + 输入指纹），使重复运行时能跳过
  未发生变化的阶段；`--force` 可绕过该机制。依赖性状的阶段按运行名分别标记。
- 运行清单（`.akgwas/run_manifest[_<run>].json`），记录工具
  版本、完整解析后的配置，以及外部工具的版本。
- 两种 multibed 后端：`python`（直接写出 `.bed`，默认）和 `plink`
  （按分块生成 `.tped` → plink）。
- 可选比对器：bowtie（`-a -v 0`）或 bwa（`mem -k31 -T31 -a -c`）。
- 对内部 K/PCA 进行确定性的分块抽样（`--seed`）。
- `--gemma-combine`（默认**关闭**）：逐个性状、全部 k-mers 的关联结果
  转储（`kmergwas_outputs/gemma_combine/`）现改为可选开启。它是单次运行中
  最大的输出，而显著 k-mers——连同其 `af`/`beta`/`se`/
  `p_wald`——已在 `sigk/` 和 `sigk_mapped/` 中，因此只有在需要重新设阈值
  或绘图时才需要完整转储。关闭时，`collect` 阶段还会跳过
  对每个分块 `.assoc.txt` 的*读取*（这本身就是一笔很大的 I/O 开销）。

### 变更

- **配置：`config.yaml` → 命令行参数（不再有配置文件）。** 每个值
  现在都是一个带默认值的参数；必填的 `--programhome` 是解析相对路径
  以及写出输出位置的基准。确切的调用方式会被记录在运行清单中。这使得
  原先临时采用的 `config.toml` 方案变得不再必要，因此 v1.0.0 发布时
  **完全不含任何配置文件**。
- **表型：每个性状一个单列表文件 → 一张多列表格。** 第一列是
  样本 ID，表头行是性状名，缺失值为 `NA`。工具会依据 bed 的 `.fam`
  校验 ID 集合，并自动重排行顺序，然后导出 GEMMA 所需的逐性状单列
  文件。这移除了旧的 `pheno_dir/` 布局及其差一（off-by-one）缺陷。
- **编排：Snakemake `.smk` + GNU parallel/rush + Bash → 纯 Python。**
  去除了 `snakemake`、`perl`、`parallel` 和 `rush` 依赖。
- **比对器默认值：`bwa mem` → `bowtie -a -v 0`**（全部全长精确命中）。
  bwa 仍可使用。
- **CLI：`--build/--kgwas/--all` 参数 → 子命令**（`akgwas build|…`）。
- **输出布局：共享 vs 按运行隔离。** `bed/`、`QK/`、`kmatrixs_out/` 与
  性状无关，是共享的；`pheno_single/`、`temp_gemma/`、
  `kmergwas_outputs/` 及其标记则按运行隔离（在设置了 `--run-name` 时位于其下）。
- 样本对齐现在是强制性的，而非建议性的：内部 K/PCA 通过
  拼接 bed 分块来构建（绝不使用 `plink --merge-list`，因为它会重排样本），
  且 PCA 样本 ID 会依据规范的 `.fam` 顺序进行检查。
- 参考索引列表：在 bowtie 默认下，`ref_list.txt` 存放的是 bowtie
  索引**前缀**（由 `bowtie-build` 生成），而非 `.sa` 路径。
- **README 的快速开始现在展示了 `bowtie-build` 步骤。** 它之前直接
  跳到 `akgwas all`，只注明 `ref_list.txt` 存放 bowtie
  索引前缀，因此首次使用的用户无从得知必须先构建索引。现在快速开始
  （以及 `ref_list.txt` 一节）展示了
  `bowtie-build <ref.fasta> <prefix>` 命令，并明确说明 akgwas 从不
  自行构建索引——由你来准备索引，并把其前缀写入
  `ref_list.txt`。

### 修复

- **`qk` 可能把 GEMMA 的*日志*当作亲缘关系矩阵保存。** `_collect_gemma_outputs`
  会在 `qk_dir/output/` 和 `qk_dir` 两处匹配 `<prefix>.*.txt`，然后取
  字典序最靠前的一个。在第二次 `qk` 运行时，上一次运行残留的
  `QK_kinship.log.txt`（位于 `qk_dir`，且 'Q' < 'o'）排在新的
  `output/QK_kinship.cXX.txt` 之前，并被移动到 `QK_kinship.txt` 之上。GEMMA 随后在
  下游拒绝该文件（"number of columns in the kinship file does not
  match the number of individuals for row = 0"），却给出误导性的
  `E602`/`.assoc.txt` 提示。现在收集器只接受真正的矩阵
  文件名——`-gk 1`（中心化，默认）对应 `<prefix>.c<digits>.txt`，
  `-gk 2`（标准化）对应 `<prefix>.s<digits>.txt`——因此 `.log.txt` 绝不会
  被误当作矩阵。该阶段会传入自身的 `-gk` 模式，以便知道预期的
  后缀。
- **当 kmtricks 运行目录已存在时，`akgwas build` 会崩溃。**
  `build` 阶段在调用 `kmtricks pipeline` 之前创建了 `kmatrixs_out/`，
  但若 `--run-dir` 已存在，kmtricks 会中止（`[CheckFailedError] ... Directory
  already exists!`）。现在该运行目录会被清空，留给 kmtricks 自行
  创建，这也让 `--force` 重跑（以及冒烟测试的强制
  第二轮）得以正常工作。
- **除非使用 `--verbose`，否则失败是静默的。** `run_module` 捕获了
  `AkgwasError` 却只写入日志文件，因此运行失败后只剩下
  一个裸的非零退出码，控制台上没有任何消息。现在 `Logger.error` 始终
  回显到 stderr，且当 `gemma -lmm` 以非零退出时，`kgwas` 阶段会追加 GEMMA
  逐性状日志的尾部内容，使数值型 LMM 失败
  （"brent did not converge"）能被看见，而不是隐藏在 `E601` 之后。
- **冒烟测试的合成数据导致 GEMMA 的 LMM 失败。** 一个 6 样本、
  2 变体的数据集会产生秩亏的 GRM：GEMMA 打印 "Matrix G has 3
  eigenvalues close to zero"，GSL 的 `brent` 优化器失败（"function value
  is not finite"），运行以非零退出。（一个完全由*稀有* k-mers 构成的
  位点集同样会失败——GRM 条件数很差。）现在 `make_synthetic`
  会植入一个宽泛的**常见** k-mers 位点集（`4 * n_samples` 个 k-mers，每个
  由 25-75% 的样本携带），且冒烟测试使用 20 个样本，因此
  GRM 条件良好，LMM 能够收敛。
- **冒烟测试的合成位点集仍会产生虚假的、与性状相关的
  k-mers。** 将每个植入的 k-mer 插入随机背景时，也会产生
  与插入片段重叠的*侧翼* k-mers。当若干携带者恰好
  共享插入位点相邻的单个碱基时，它们也就共享一个侧翼
  k-mer；该侧翼变异是植入 k-mer 携带者的严格子集，因此它与性状相关，
  并同时污染关联检验和 GRM。失败的那次运行显示出
  `af=0.10` 的显著 k-mers（即 2/20 个携带者），这是位点集设计
  绝不可能产生的，而植入的因果 k-mer 反而落出了 `sigK`。现在 `make_synthetic`
  使每个**因果** k-mer 最多有四个携带者，并为每个携带者写出**互不相同**的
  侧翼碱基，因此每个侧翼 k-mer 都唯一对应一个基因组，并被
  `recurrence-min` 丢弃（它只出现一次）。只有因果 k-mers 需要这样处理：
  常见位点集 k-mers 的侧翼与性状无关，因此无害。每个 k-mer 在植入前
  也会被规范化，因此无论计数器如何处理链方向，运行都是
  可复现的。
- 冒烟测试现在接受任意一条链上的植入 k-mer（k-mer 计数器
  可能存储规范化后的形式）。
- 表型/fof/样本顺序不匹配曾静默地产生错误的关联结果：现在会
  校验 ID 集合，并将行重排为 `.fam` 顺序。
- 构建显著 k-mers 的 FASTA 时，第一个显著 k-mer 可能因 `NR>1`
  差一而被丢弃——现在 k-mer 直接取自
  sigK 行。
- `--multibed-threads` 被忽略：流式的 `kmtricks aggregate` 步骤现在会
  遵循它（此前用的是 `--build-threads`）；`--qk-threads` 现在会
  传给 `plink --pca`（此前被忽略）；已废弃的 `--kgwas-threads`
  参数被删除（GEMMA 作业是单线程的——请使用 `--cores`）。
- 当 `kmtricks aggregate` 向其 stderr 写入超过一个管道缓冲区的数据、
  而其 stdout 同时正被消费时，可能发生死锁——现在 stderr 会写入
  临时文件。
- `readme` 文件名大小写（`Readme.md` → `README.md`），以便构建后端在
  区分大小写的文件系统上能找到它。
- `recipe/meta.yaml` 的 `about.home` 指向了错误的仓库名。

### 移除

- 配置文件完全移除（`config.yaml` / 任何临时的 `config.toml`）：
  所有设置现在都是命令行参数。
- 运行时依赖 `snakemake`、`perl`、`parallel`、`rush`。
- 中间的多 GB 文本矩阵和 `.tped` 文件：矩阵由
  `kmtricks aggregate` 流式输出，且在默认后端下，绝不写出任何
  `.tped`。

## [0.1.0]

初始版本：一个由 Bash 脚本驱动的 Snakemake 工作流，驱动 kmtricks、plink、GEMMA
和 bwa/samtools。配置位于 `config.yaml`；每个性状一个
单列表型文件，放在 `pheno_dir/` 下。
