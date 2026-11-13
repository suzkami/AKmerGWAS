# 验证清单（在 Linux 机器上运行）

此代码树是在 Windows 上编写的，当时无法执行它。下面所有内容都应当在将要
构建/测试 conda 包的那台机器上运行**一次**。请自上而下依次操作；每一步都
足够独立，因此某一步失败就能定位到具体的层。

## 0. 环境

```bash
conda create -n akgwas -c conda-forge -c bioconda \
    python=3.12 kmtricks gemma plink bowtie bwa pytest
conda activate akgwas
cd v6
```

## 1. 静态检查（无需外部工具）

```bash
python -m compileall akgwas tests          # every module must byte-compile
python -c "import akgwas, akgwas.stages"   # import graph is clean
```

然后是单元测试：

```bash
pytest tests/test_units.py -q
```

这会检验纯逻辑部分：大小解析、bed 打包/解包往返、bed 拼接/取子集、
表型重排/导出 + ID 不匹配、MAF 过滤 + 分块、SAM 过滤、fof 解析、
配置默认值 + 校验。

要特别留意：

- `test_bed_concat_and_subset`——确认我们手工实现的 `.bed` 字节布局
  与 PLINK 读取的一致。如果它失败，说明 `python` 后端有误。
- `test_pheno_reorder_and_export`——确认 ID→fam 的重排。

## 2. CLI 界面

```bash
python -m akgwas --version          # -> akgwas 1.0.0
python -m akgwas --errors           # prints the E1xx..E9xx table
python -m akgwas build --help       # grouped flags render
python -m akgwas all --help
python -m akgwas all                # must fail: --programhome is required
```

## 3. 合成数据端到端

```bash
python -m tests.make_synthetic /tmp/syn
# builds a bowtie index and runs the whole pipeline, then re-runs with --force
pytest tests/test_smoke.py -q
```

冒烟测试断言：植入的 `traitA` k-mer 出现在
`sigk/traitA_sigK.txt` 中，`siggeno`/`samfile`/`sigk_mapped` 输出存在，并且
强制重跑的结果是**字节完全一致**的（可重复性）。如果缺少 bowtie，
该测试会被跳过——用 (`-rs`) 检查它是否真的运行了。

这里落实了三项修复。数据集在 20 个样本上植入了一个宽泛的**常见** k-mer
面板（不是 6 个样本 / 2 个变异，也不是一个稀有 k-mer 面板），
这样 GEMMA 的 LMM 能够收敛，而不会在秩亏或条件数很差的 GRM 上于 GSL `brent`
中崩溃。只有两个**因果** k-mer 保持稀疏（至多四个携带者）并带有
**不同的侧翼碱基**，因此因果 k-mer 的任何侧翼 k-mer 都无法通过
`recurrence-min` 存留下来与性状相关（此前这类侧翼既污染了 GRM，又把真正的
因果 k-mer 挤出了 `sigK`）；常见面板 k-mer 的侧翼与性状无关，是无害的。
此外，强制重跑会在 `kmatrixs_out/` 已存在的情况下重新进入 `build` 阶段——
该路径正是本次修订修复的 kmtricks 运行目录 bug（该阶段现在会先清空，
再让 kmtricks 重新创建该目录）。

如果运行失败，原因会被打印到 stderr（对于 `E601` 还会附带 GEMMA 日志
末尾），而不会被吞掉。

## 4. 打包

```bash
# from the repo root that will be tagged
python -m pip install build && python -m build      # sdist + wheel
python -m pip install dist/akgwas-1.0.0-py3-none-any.whl
akgwas --version
```

Conda recipe 部分：

```bash
# after tagging v1.0.0 and filling recipe/meta.yaml's sha256:
conda-build v6/recipe/meta.yaml -c conda-forge -c bioconda
```

在进行上述操作时需要确认的注意事项：

- `readme = "README.md"`（大写 R）能够正确解析——0.1.0 的构建曾因
  小写的 `Readme.md` 而失败。
- wheel 中不包含任何数据文件——该包是纯 Python，运行时不会从 `akgwas/`
  读取任何内容（`packages=` 中没有配置模板）。
- recipe 会拉取 **bowtie**（新的默认比对器），并去掉
  snakemake/perl/rush/parallel。

## 5. 只有人工才能完成的接口检查

- **外部 K/PCA 路径。** 传入 `--use-external-qk --kinship K --pca PCA`，
  并提供一个你信任的真实亲缘关系 + PCA，确认重排行为符合预期。分别尝试
  不带 `--pca-has-id`（假定按行顺序）和带 `--pca-has-id`（对 ID
  进行重排/校验）两种情况。
- **bwa 路径。** 加上 `--aligner bwa`，让 `ref_list.txt` 指向 bwa 索引前缀，
  仅重新运行 `akgwas alignment --programhome ... --force`。
- **plink multibed 后端。** 加上 `--bed-backend plink`，确认各分块与
  `python` 后端完全一致（对两者分别执行 `plink --bfile matrix_0001 --recode`，
  再做 diff）。
- **真实耐受性。** 在一个小型真实数据集上运行，检查显著 k-mer 是否
  映射到合理的参考（基因组）位置。

## 已知的待审点（请仔细审查）

1. **GEMMA `.assoc.txt` 的列顺序**——`kgwas.py` 假定为标准 12
   列，并取用索引 `[1,6,7,8,11]`（rs、af、beta、se、p_wald）。请对照
   已安装的 GEMMA 版本加以确认。
2. **gemma `-gk` 的输出文件名**——`qk.py` 现在只接受真正的矩阵文件名：
   `-gk 1`（中心化，默认值）对应 `QK_kinship.cXX.txt`，`-gk 2`（标准化）
   对应 `QK_kinship.sXX.txt`。它绝不会误取 `QK_kinship.log.txt`。请确认
   你的 GEMMA 构建使用的正是这些后缀。
3. **bowtie 全长精确匹配过滤器**——`-a -v 0 -f -S`，外加 `alignment.py` 中
   对 CIGAR `== kM` 和 `NM:i:0` 的检查。请确认你所用版本的 bowtie 对精确
   命中会输出 `NM:i:0`。
4. **plink `--pca` 特征向量的形状**——`qk.py` 期望 FID IID PC1..PC3 且无
   表头。请确认你的 plink 版本。
