# 故障排查

流水线能够归因的每一种失败都会打印一个稳定的错误码，例如

```
[E301] phenotype sample IDs do not match the bed .fam sample set
    hint: Phenotype/fof sample IDs must form the same set as the bed .fam IDs.
```

随时可以用 `akgwas --errors` 打印完整的错误码表。错误码按阶段分组：

| 范围 | 领域 |
|---|---|
| `E1xx` | 环境 / 预检（工具、ulimit、路径） |
| `E2xx` | 输入校验（fof、表型、参考列表、配置） |
| `E3xx` | 样本对齐 / 样本顺序不变量 |
| `E4xx` | build + multibed（kmtricks、分块、bed 写入） |
| `E5xx` | 亲缘关系（kinship）/ PCA（plink、gemma -gk） |
| `E6xx` | GWAS（gemma -lmm） |
| `E7xx` | 收集 / 合并 |
| `E8xx` | 比对（bowtie / bwa / SAM 解析） |
| `E9xx` | 内部 / 未预期 |

---

## 错误码参考

### E1xx — 环境 / 预检

| 错误码 | 含义 | 修复方法 |
|---|---|---|
| `E101` | 某个必需的外部工具不在 `PATH` 上 | 激活 conda 环境，或执行 `conda install -c bioconda <tool>` |
| `E102` | 打开文件数上限（`ulimit -n`）低于 `--ulimit-min` | 运行前执行 `ulimit -n 4096`，或调低 `--ulimit-min` |
| `E103` | 工作/输出目录不可写 | 检查权限和磁盘剩余空间 |

### E2xx — 输入校验

| 错误码 | 含义 | 修复方法 |
|---|---|---|
| `E201` | 某个输入文件不存在 | 检查命令行参数（flag）/路径；相对路径以 `--programhome` 为基准解析 |
| `E202` | 无效设置（取值错误、块大小范围错误、未知选项） | 检查该取值；`akgwas <subcommand> --help` 会列出有效范围 |
| `E204` | `fof.txt` 格式错误 | 每一行必须是 `name : file1 ; file2 ; ...` |
| `E205` | 行数不匹配（表型 vs fof） | 每个 fof 样本对应一行表型 |
| `E206` | 参考列表为空 | 向 `ref_list.txt` 中至少添加一行索引 |
| `E207` | 表型表格式错误（各行长度不一致、缺少表头） | 第一行为表头，第一列为样本 ID |

### E3xx — 样本对齐（不变量）

| 错误码 | 含义 | 修复方法 |
|---|---|---|
| `E301` | 表型/fof 的样本 ID 集合 ≠ bed `.fam` 的 ID 集合 | 对齐样本集合；检查拼写错误、空白字符、重复 ID |
| `E302` | K/PCA 的样本顺序 ≠ GWAS bed 的顺序 | 在相同的样本顺序上重建 K/PCA（见下文） |
| `E303` | 无法校验样本顺序（PCA 不含 ID） | 如果 PCA 文件带有 ID，请加上 `--pca-has-id` |

### E4xx — build + multibed

| 错误码 | 含义 | 修复方法 |
|---|---|---|
| `E401` | `kmtricks` 的 count/merge 失败 | 查看 `.akgwas/` 中捕获的 kmtricks 日志 |
| `E402` | kmtricks 生成了空矩阵 | 输入为空，或 `--hard-min`/`--recurrence-min` 过于严格 |
| `E403` | 未知的 bed 后端 | 使用 `--bed-backend python` 或 `--bed-backend plink` |
| `E404` | 矩阵某行格式错误（宽度不对 / 非 0-1） | 通常是中间输入损坏；用 `--force` 重新运行 `build` |

### E5xx — 亲缘关系（kinship）/ PCA

| 错误码 | 含义 | 修复方法 |
|---|---|---|
| `E501` | `plink --pca` 失败 | 检查 plink 日志；通常是 bed 格式错误 |
| `E502` | `gemma -gk` 失败 | 检查 gemma 日志；通常是亲缘关系输入的问题 |

### E6xx — GWAS

| 错误码 | 含义 | 修复方法 |
|---|---|---|
| `E601` | 某个分块×性状的 `gemma -lmm` 失败 | 查看 `logs/` 下按性状分文件的 gemma 日志（日志末尾也会附加到错误信息中） |
| `E602` | 某个 GEMMA `.assoc.txt` 的列数不符合预期 | 检查 GEMMA 版本；预期为 12 列 |

> **数值型 LMM 失败**（"Matrix G has N eigenvalues close to zero"、"GSL
> ERROR: ... brent.c"、"Brent did not converge"）会以 `E601` 报告。它们
> 意味着亲缘关系矩阵秩亏——几乎总是样本太少和/或 k-mer 变异太少
>（例如一个很小的试点数据集）。扩大面板（增加样本、调低
> `--hard-min`/`--recurrence-min`、调低 `--maf`），或通过 `--use-external-qk`
> 提供一个外部良态（well-conditioned）的亲缘关系。

### E7xx — 收集

| 错误码 | 含义 | 修复方法 |
|---|---|---|
| `E701` | 没有可收集的内容（没有分块，或没有显著的 k-mer） | 若无命中，这只是提示信息；调低 `--threshold` |

### E8xx — 比对

| 错误码 | 含义 | 修复方法 |
|---|---|---|
| `E801` | 未知的比对器 | 使用 `--aligner bowtie` 或 `--aligner bwa` |
| `E802` | 比对器失败 | 检查比对器的 stderr；检查 `ref_list.txt` 中的索引前缀 |

### E9xx — 内部错误

| 错误码 | 含义 | 修复方法 |
|---|---|---|
| `E999` | 未预期的内部错误 | 连同 `.akgwas/run_manifest.json` 一起报告 |

---

## 样本顺序不变量（如果你遇到 E301/E302/E303，请阅读本节）

一切都以 **bed `.fam` 的顺序**为准，而该顺序来自
`fof.txt`（样本声明顺序——*并非*排序后的顺序）。三处输入必须一致：

1. **表型**——按 ID 与 `.fam` 匹配，然后自动重排。如果 ID *集合*
   不同，就会得到 `E301`。
2. **亲缘关系（kinship）/ PCA**——必须在与最终 GWAS bed *相同*的样本顺序
   上构建。不使用 `--use-external-qk` 时，工具会自行拼接 bed 分块
   （绝不使用 `plink --merge-list`，它会按字典序静默重排样本）来构建，
   并校验 PCA 的 ID，因此在工具内部你不可能弄错。
3. **外部 K/PCA**——使用 `--use-external-qk` 时，顺序由*你*负责。
   - 如果你的 PCA 文件在前几列带有样本 ID，请加上
     `--pca-has-id`；工具会据此重排并校验。
   - 亲缘关系矩阵**不含 ID**，因此假定它遵循你提供给 PCA 的相同顺序。
     如有疑问，请从同一 SNP 集合、以相同的样本顺序重新一起生成 K 和 PCA。

GWAS bed 的 `.fam` 是唯一可信来源。要打印它：

```bash
cut -d' ' -f2 bed/matrix_0001.fam
```

将那一列与你的外部 PCA ID 进行比对——它们必须完全一致。

---

## 可重复性与续跑

- 每个阶段都会在 `.akgwas/` 下写入一个原子化的 `<stage>.done` 完成标记
  （附带输入指纹）。重新运行时，输入未变的阶段会被跳过。
- `--force` 会忽略这些标记，重新运行所有内容。
- 依赖性状的阶段（`kgwas`、`collect`、`alignment`）按 `--run-name` 分别
  标记（`<run-name>.<stage>.done`），因此共享同一个 `outdir` 的不同表型运行
  既不会互相覆盖，也不会错误地互相跳过。
- 抽样分块的随机种子是 `--seed`；工具版本、完整解析后的配置以及所选工具的
  版本都记录在 `.akgwas/run_manifest[_<run-name>].json` 中。请将
  该文件与你的结果一并保存。

## 没有任何输出 / 试运行

`--dry-run` 会打印它*将要*执行的每一条外部命令，但不做任何改动。
用它可以在耗时的构建之前检查路径是否合理。
