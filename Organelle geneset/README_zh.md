# 人类细胞器基因集

收集自多个权威数据库的人类细胞器基因集，经整理后每个基因集仅保留各细胞器的独有基因。

## 数据来源

| 细胞器 | 数据库 | 原始基因数 |
|--------|--------|-----------|
| 线粒体 Mitochondria | [MitoCarta 3.0](https://www.broadinstitute.org/mitocarta) | 1077 |
| 内质网 Endoplasmic Reticulum | [MsigDB C5:GO:CC](https://www.gsea-msigdb.org/gsea/msigdb) | 1652 |
| 高尔基体 Golgi | [MsigDB C5:GO:CC](https://www.gsea-msigdb.org/gsea/msigdb) | 1586 |
| 溶酶体 Lysosome | [hLGDB](http://lysosome.unipg.it/) | 435 |
| 过氧化物酶体 Peroxisome | [PeroxisomeDB](https://www.peroxisomedb.org/) | 98 |
| 脂滴 Lipid Droplet | [Lipid Droplet Knowledge Portal](https://lipiddroplet.org/) | 1077 (≥2) / 505 (≥3) |

## 数据处理

1. **跨细胞器重复基因去除**：在多个细胞器基因集中同时出现的基因已被删除，每个基因集仅保留该细胞器的独有基因。
2. **脂滴基因集筛选**：设置了两个置信度阈值——在 ≥2 个研究中检测到的基因（min2，550个）和更严格的 ≥3 个研究阈值（min3，232个）。下游处理分析使用 min3 版本。

## 基因数量（处理后）

| 编号 | 基因集 | 独有基因数 |
|------|--------|-----------|
| 00 | 全部细胞器 All organelles | 3934 |
| 01 | 内质网 Endoplasmic Reticulum | 903 |
| 02 | 高尔基体 Golgi | 975 |
| 03 | 脂滴 Lipid Droplet (min2) | 550 |
| 04 | 脂滴 Lipid Droplet (min3) | 232 |
| 05 | 溶酶体 Lysosome | 238 |
| 06 | 线粒体 Mitochondria | 988 |
| 07 | 过氧化物酶体 Peroxisome | 48 |
