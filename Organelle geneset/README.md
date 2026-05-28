# Human Organelle Gene Sets

Curated human gene sets for 6 subcellular organelles, collected from established databases and processed to retain only organelle-unique genes.

## Sources

| Organelle | Database | Original Genes |
|-----------|----------|---------------|
| Mitochondria | [MitoCarta 3.0](https://www.broadinstitute.org/mitocarta) | 1077 |
| Endoplasmic Reticulum | [MsigDB C5:GO:CC](https://www.gsea-msigdb.org/gsea/msigdb) | — |
| Golgi | [MsigDB C5:GO:CC](https://www.gsea-msigdb.org/gsea/msigdb) | — |
| Lysosome | [hLGDB](http://lysosome.unipg.it/) | — |
| Peroxisome | [PeroxisomeDB](https://www.peroxisomedb.org/) | — |
| Lipid Droplet | [Lipid Droplet Knowledge Portal](https://lipiddroplet.org/) | 1077 (≥2) / 505 (≥3) |

## Processing

1. **Cross-organelle overlap removal**: Genes shared across multiple organelle gene sets were removed. Each gene set retains only organelle-unique genes.
2. **Lipid droplet filtering**: Two confidence thresholds were applied — genes detected in ≥2 studies (min2, 550 genes) and a more stringent ≥3 studies threshold (min3, 232 genes). The min3 set is recommended for downstream analysis.

## Gene Counts (after processing)

| # | Gene Set | Unique Genes |
|---|----------|-------------|
| 00 | All organelles | 3934 |
| 01 | Endoplasmic Reticulum | 903 |
| 02 | Golgi | 975 |
| 03 | Lipid Droplet (min2) | 550 |
| 04 | Lipid Droplet (min3) | 232 |
| 05 | Lysosome | 238 |
| 06 | Mitochondria | 988 |
| 07 | Peroxisome | 48 |
