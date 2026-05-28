"""

Dependencies
------------
    pip install numpy pandas scikit-learn scipy
    # optional but recommended:
    pip install pyarrow shap matplotlib statsmodels

Author note: v1.3 beta version
"""

from __future__ import annotations

import os
import json
import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from sklearn.base import clone
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=UserWarning)
RNG = np.random.default_rng(42)


# ======================================================================================
# 1. CONFIGURATION  
# ======================================================================================
@dataclass
class Config:
    # ---- inputs ----
    organelle_geneset_path: str = "Organelle geneset"   # GMT file 
    gene_matrix_path: str = "data/ukb_gene_level_matrix.parquet"  # rows=eid, cols=gene symbols
    phenotype_path: str = "data/ukb_phenotypes.parquet"           # eid + covariates + raw UKB fields
    eid_col: str = "eid"

    # ---- covariates (adjust to your extract's column names) ----
    covariate_cols: Sequence[str] = field(default_factory=lambda: [
        "age", "sex", "bmi",
        "pc1", "pc2", "pc3", "pc4", "pc5",
        "pc6", "pc7", "pc8", "pc9", "pc10",
    ])

    # ---- outcome to model (pick one key from DISEASE_ICD below, or "any_metabolic") ----
    outcome: str = "t2d"

    # ---- UKB field names holding diagnosis codes  ----
    icd10_cols_prefix: str = "icd10"   # any columns starting with this are treated as ICD-10 arrays
    icd9_cols_prefix: str = "icd9"
    self_report_cols_prefix: str = "selfreport_noncancer"  # UKB field 20002 array

    # ---- organelle-score aggregation method ----
    score_method: str = "zmean"   # one of: "zmean", "rank", "pca1"
    min_genes_per_organelle: int = 3      # drop organelle sets with too few mapped genes
    standardize_genes: bool = True

    # ---- interaction features ----
    include_pairwise_interactions: bool = True
    include_covariates_in_model: bool = True

    # ---- modelling ----
    n_splits: int = 5
    test_fraction: float = 0.2
    n_bootstrap: int = 1000
    class_weight: Optional[str] = "balanced"

    # ---- outputs ----
    output_dir: str = "outputs"


# Algorithmically-defined metabolic outcomes via ICD code prefixes (extend as needed).
# Self-report (UKB 20002) codes added where commonly used.
DISEASE_ICD: Dict[str, Dict[str, List[str]]] = {
    "t2d":          {"icd10": ["E11"],                 "icd9": ["250"], "selfreport": ["1223"]},
    "dyslipidemia": {"icd10": ["E78"],                 "icd9": ["272"], "selfreport": ["1473"]},
    "hypertension": {"icd10": ["I10", "I15"],          "icd9": ["401"], "selfreport": ["1065", "1072"]},
    "ihd":          {"icd10": ["I20", "I21", "I22",
                               "I23", "I24", "I25"],   "icd9": ["410", "411", "412", "413", "414"],
                                                       "selfreport": ["1075"]},
    "masld":        {"icd10": ["K760", "K758"],        "icd9": ["5715"], "selfreport": []},
    "obesity":      {"icd10": ["E66"],                 "icd9": ["278"], "selfreport": ["1073"]},
}


# ======================================================================================
# 2. LOAD ORGANELLE GENE SETS
# ======================================================================================
def load_organelle_genesets(path: str) -> Dict[str, set]:
    """
    Load organelle -> set(gene symbols).

    Supports:
      * a single .gmt file  (name<TAB>desc<TAB>gene1<TAB>gene2 ...)
      * a directory of per-organelle files (.txt/.csv/.tsv); filename = organelle name,
        one gene per line (or first column).
    """
    genesets: Dict[str, set] = {}

    if os.path.isfile(path) and path.lower().endswith(".gmt"):
        with open(path, "r") as fh:
            for line in fh:
                parts = [p.strip() for p in line.rstrip("\n").split("\t") if p.strip()]
                if len(parts) >= 3:
                    name, _desc, *genes = parts
                    genesets[name] = {g.upper() for g in genes}
    elif os.path.isdir(path):
        for fname in sorted(os.listdir(path)):
            fpath = os.path.join(path, fname)
            if not os.path.isfile(fpath):
                continue
            organelle = os.path.splitext(fname)[0]
            # Read line-by-line and take the first token on each line. We split only on
            # comma / tab (never let pandas "sniff" a delimiter, which can shatter symbols
            # like GENE0 on the letter E). Real HGNC symbols contain no comma/tab.
            tokens: List[str] = []
            with open(fpath) as fh:
                for ln in fh:
                    ln = ln.strip()
                    if not ln:
                        continue
                    first = ln.replace("\t", ",").split(",")[0].strip().strip('"').strip("'")
                    if first:
                        tokens.append(first)
            genes = {t.upper() for t in tokens if t.lower() not in {"gene", "gene_symbol", "symbol"}}
            if genes:
                genesets[organelle] = genes
    else:
        raise FileNotFoundError(f"Organelle gene-set path not found: {path}")

    if not genesets:
        raise ValueError(f"No gene sets parsed from {path}")
    print(f"[genesets] loaded {len(genesets)} organelle sets: {list(genesets)}")
    return genesets


# ======================================================================================
# 3. LOAD GENE-LEVEL MATRIX + PHENOTYPES
# ======================================================================================
def _read_table(path: str) -> pd.DataFrame:
    if path.lower().endswith(".parquet"):
        return pd.read_parquet(path)
    if path.lower().endswith((".tsv", ".txt")):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)


def load_gene_matrix(cfg: Config) -> pd.DataFrame:
    """rows = individuals (indexed by eid), columns = gene symbols (float values)."""
    df = _read_table(cfg.gene_matrix_path)
    if cfg.eid_col in df.columns:
        df = df.set_index(cfg.eid_col)
    df.columns = [str(c).upper() for c in df.columns]
    df = df.apply(pd.to_numeric, errors="coerce")
    print(f"[gene matrix] {df.shape[0]} individuals x {df.shape[1]} genes")
    return df


def load_phenotypes(cfg: Config) -> pd.DataFrame:
    df = _read_table(cfg.phenotype_path)
    if cfg.eid_col in df.columns:
        df = df.set_index(cfg.eid_col)
    print(f"[phenotypes] {df.shape[0]} individuals x {df.shape[1]} fields")
    return df


# ======================================================================================
# 4. BUILD BINARY DISEASE OUTCOMES FROM UKB DIAGNOSIS FIELDS
# ======================================================================================
def _match_any_prefix(values: pd.Series, prefixes: Sequence[str]) -> pd.Series:
    """True where the concatenated code string starts with any given prefix."""
    if not prefixes:
        return pd.Series(False, index=values.index)
    s = values.fillna("").astype(str).str.upper().str.replace(".", "", regex=False)
    mask = pd.Series(False, index=values.index)
    for p in prefixes:
        mask |= s.str.startswith(p.upper())
    return mask


def build_outcome(phen: pd.DataFrame, cfg: Config, outcome_key: str) -> pd.Series:
    """
    Construct a binary case/control label for one disease by scanning all ICD-10/ICD-9/
    self-report array columns for the configured code prefixes.
    """
    if outcome_key == "any_metabolic":
        sub = pd.DataFrame(index=phen.index)
        for k in DISEASE_ICD:
            sub[k] = build_outcome(phen, cfg, k)
        return (sub.sum(axis=1) > 0).astype(int).rename("any_metabolic")

    spec = DISEASE_ICD[outcome_key]
    icd10_cols = [c for c in phen.columns if str(c).lower().startswith(cfg.icd10_cols_prefix)]
    icd9_cols = [c for c in phen.columns if str(c).lower().startswith(cfg.icd9_cols_prefix)]
    sr_cols = [c for c in phen.columns if str(c).lower().startswith(cfg.self_report_cols_prefix)]

    label = pd.Series(False, index=phen.index)
    for c in icd10_cols:
        label |= _match_any_prefix(phen[c], spec.get("icd10", []))
    for c in icd9_cols:
        label |= _match_any_prefix(phen[c], spec.get("icd9", []))
    for c in sr_cols:
        label |= _match_any_prefix(phen[c], spec.get("selfreport", []))

    y = label.astype(int).rename(outcome_key)
    print(f"[outcome:{outcome_key}] cases={int(y.sum())} ({100*y.mean():.2f}%)  N={len(y)}")
    return y


# ======================================================================================
# 5. PER-ORGANELLE SCORES
# ======================================================================================
def compute_organelle_scores(
    gene_matrix: pd.DataFrame,
    genesets: Dict[str, set],
    cfg: Config,
) -> pd.DataFrame:
    """
    Aggregate gene-level values into one score per organelle per individual.

    Methods:
      zmean : z-score each gene, then mean across the organelle's mapped genes
      rank  : per-individual rank-normalise genes, then mean rank (ssGSEA-like, lightweight)
      pca1  : first principal component of the organelle's genes (sign-aligned to mean)
    """
    X = gene_matrix.copy()

    if cfg.score_method == "rank":
        # rank within each individual across all genes, scaled to [0,1]
        X = X.rank(axis=1, pct=True)
    elif cfg.standardize_genes:
        mu = X.mean(axis=0)
        sd = X.std(axis=0).replace(0, np.nan)
        X = (X - mu) / sd

    scores = {}
    for organelle, genes in genesets.items():
        cols = [g for g in genes if g in X.columns]
        if len(cols) < cfg.min_genes_per_organelle:
            print(f"  [skip] {organelle}: only {len(cols)} genes mapped "
                  f"(< {cfg.min_genes_per_organelle})")
            continue
        block = X[cols]
        if cfg.score_method == "pca1":
            from sklearn.decomposition import PCA
            filled = block.fillna(block.mean())
            comp = PCA(n_components=1, random_state=42).fit_transform(filled).ravel()
            # align sign so higher score ~ higher mean expression/burden
            if np.corrcoef(comp, filled.mean(axis=1))[0, 1] < 0:
                comp = -comp
            scores[organelle] = comp
        else:
            scores[organelle] = block.mean(axis=1, skipna=True).values

    out = pd.DataFrame(scores, index=X.index)
    out = out.apply(lambda c: (c - c.mean()) / c.std() if c.std() else c)  # standardise scores
    print(f"[organelle scores] {out.shape[1]} organelle scores for {out.shape[0]} individuals")
    return out


# ======================================================================================
# 6. INTERACTION FEATURES  (the core "organelle interaction" object)
# ======================================================================================
def build_interaction_features(org_scores: pd.DataFrame, cfg: Config) -> pd.DataFrame:
    """Main organelle effects + (optionally) all pairwise organelle x organelle products."""
    feats = org_scores.copy()
    if cfg.include_pairwise_interactions:
        cols = list(org_scores.columns)
        for i in range(len(cols)):
            for j in range(i + 1, len(cols)):
                a, b = cols[i], cols[j]
                feats[f"INT__{a}__x__{b}"] = org_scores[a] * org_scores[b]
    print(f"[features] {feats.shape[1]} features "
          f"({org_scores.shape[1]} main + {feats.shape[1] - org_scores.shape[1]} interactions)")
    return feats


def assemble_design_matrix(
    org_feats: pd.DataFrame,
    phen: pd.DataFrame,
    y: pd.Series,
    cfg: Config,
) -> Tuple[pd.DataFrame, pd.Series, List[str]]:
    """Join organelle features + covariates, align on eid, drop rows with no outcome."""
    parts = [org_feats]
    cov_cols: List[str] = []
    if cfg.include_covariates_in_model:
        present = [c for c in cfg.covariate_cols if c in phen.columns]
        missing = [c for c in cfg.covariate_cols if c not in phen.columns]
        if missing:
            print(f"  [warn] covariates not found, skipped: {missing}")
        cov = phen[present].copy()
        if "sex" in cov.columns and cov["sex"].dtype == object:
            cov["sex"] = (cov["sex"].astype(str).str.lower()
                          .map({"male": 1, "m": 1, "1": 1, "female": 0, "f": 0, "0": 0}))
        cov_cols = present
        parts.append(cov)

    X = pd.concat(parts, axis=1, join="inner")
    common = X.index.intersection(y.index)
    X, y = X.loc[common], y.loc[common]
    X = X.replace([np.inf, -np.inf], np.nan)
    print(f"[design] {X.shape[0]} individuals x {X.shape[1]} features; "
          f"prevalence={100*y.mean():.2f}%")
    return X, y, cov_cols


# ======================================================================================
# 7. MODELS
# ======================================================================================
def make_models(cfg: Config) -> Dict[str, Pipeline]:
    """
    Three complementary models:
      logit_interaction : L2 logistic on organelle main + interaction terms (interpretable;
                          interaction coefficients are the scientific deliverable)
      logit_elasticnet  : sparse elastic-net logistic (selects informative interactions)
      gbm               : HistGradientBoosting (captures nonlinear organelle interplay)
    """
    pre = Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
    ])
    models = {
        "logit_interaction": Pipeline([
            ("pre", pre),
            ("clf", LogisticRegression(penalty="l2", C=1.0, max_iter=5000,
                                       class_weight=cfg.class_weight)),
        ]),
        "logit_elasticnet": Pipeline([
            ("pre", clone(pre)),
            ("clf", LogisticRegression(penalty="elasticnet", solver="saga",
                                       l1_ratio=0.5, C=0.5, max_iter=8000,
                                       class_weight=cfg.class_weight)),
        ]),
        "gbm": Pipeline([
            ("impute", SimpleImputer(strategy="median")),
            ("clf", HistGradientBoostingClassifier(
                learning_rate=0.05, max_iter=400, max_leaf_nodes=31,
                l2_regularization=1.0, early_stopping=True, random_state=42)),
        ]),
    }
    return models


# ======================================================================================
# 8. EVALUATION
# ======================================================================================
def _bootstrap_ci(y_true: np.ndarray, y_prob: np.ndarray, metric, n: int) -> Tuple[float, float, float]:
    point = metric(y_true, y_prob)
    idx = np.arange(len(y_true))
    vals = []
    for _ in range(n):
        s = RNG.choice(idx, size=len(idx), replace=True)
        if len(np.unique(y_true[s])) < 2:
            continue
        vals.append(metric(y_true[s], y_prob[s]))
    lo, hi = np.percentile(vals, [2.5, 97.5]) if vals else (np.nan, np.nan)
    return float(point), float(lo), float(hi)


def evaluate_cv(model: Pipeline, X: pd.DataFrame, y: pd.Series, cfg: Config) -> Dict:
    cv = StratifiedKFold(n_splits=cfg.n_splits, shuffle=True, random_state=42)
    prob = cross_val_predict(model, X.values, y.values, cv=cv,
                             method="predict_proba", n_jobs=-1)[:, 1]
    yt = y.values
    auroc = _bootstrap_ci(yt, prob, roc_auc_score, cfg.n_bootstrap)
    auprc = _bootstrap_ci(yt, prob, average_precision_score, cfg.n_bootstrap)
    brier = brier_score_loss(yt, prob)
    return {"auroc": auroc, "auprc": auprc, "brier": brier, "oof_prob": prob}


# ======================================================================================
# 9. INTERPRETATION
# ======================================================================================
def interpret_interactions(fitted: Pipeline, feature_names: Sequence[str], top_k: int = 25) -> pd.DataFrame:
    """Extract standardised coefficients; surface the strongest organelle-organelle terms."""
    clf = fitted.named_steps.get("clf")
    if not hasattr(clf, "coef_"):
        return pd.DataFrame()
    coefs = clf.coef_.ravel()
    df = pd.DataFrame({"feature": list(feature_names), "coef": coefs})
    df["odds_ratio"] = np.exp(df["coef"])
    df["abs"] = df["coef"].abs()
    df["is_interaction"] = df["feature"].str.startswith("INT__")
    return df.sort_values("abs", ascending=False).head(top_k).reset_index(drop=True)


def shap_summary(fitted_gbm: Pipeline, X: pd.DataFrame, out_png: str) -> bool:
    try:
        import shap
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:  # optional dependency
        print(f"  [shap] skipped ({e})")
        return False
    Xt = fitted_gbm.named_steps["impute"].transform(X.values)
    expl = shap.TreeExplainer(fitted_gbm.named_steps["clf"])
    sv = expl.shap_values(Xt)
    shap.summary_plot(sv, Xt, feature_names=list(X.columns), show=False, max_display=20)
    plt.tight_layout(); plt.savefig(out_png, dpi=150, bbox_inches="tight"); plt.close()
    print(f"  [shap] saved -> {out_png}")
    return True


# ======================================================================================
# 10. ORCHESTRATION
# ======================================================================================
def run(cfg: Config) -> Dict:
    os.makedirs(cfg.output_dir, exist_ok=True)

    genesets = load_organelle_genesets(cfg.organelle_geneset_path)
    gene_matrix = load_gene_matrix(cfg)
    phen = load_phenotypes(cfg)

    y = build_outcome(phen, cfg, cfg.outcome)
    org_scores = compute_organelle_scores(gene_matrix, genesets, cfg)
    org_feats = build_interaction_features(org_scores, cfg)
    X, y, _cov_cols = assemble_design_matrix(org_feats, phen, y, cfg)

    models = make_models(cfg)
    results = {}
    for name, model in models.items():
        print(f"\n=== {name} ===")
        res = evaluate_cv(model, X, y, cfg)
        a = res["auroc"]; p = res["auprc"]
        print(f"  AUROC {a[0]:.3f} [{a[1]:.3f}-{a[2]:.3f}]   "
              f"AUPRC {p[0]:.3f} [{p[1]:.3f}-{p[2]:.3f}]   Brier {res['brier']:.4f}")
        results[name] = {k: v for k, v in res.items() if k != "oof_prob"}

    # Fit interaction model on full data and report organelle-organelle effects
    fitted = models["logit_interaction"].fit(X.values, y.values)
    inter = interpret_interactions(fitted, X.columns, top_k=30)
    if not inter.empty:
        inter.to_csv(os.path.join(cfg.output_dir, f"{cfg.outcome}_organelle_interaction_effects.csv"),
                     index=False)
        print("\nTop organelle / interaction effects:")
        print(inter[["feature", "coef", "odds_ratio", "is_interaction"]].head(15).to_string(index=False))

    fitted_gbm = models["gbm"].fit(X.values, y.values)
    shap_summary(fitted_gbm, X, os.path.join(cfg.output_dir, f"{cfg.outcome}_gbm_shap.png"))

    with open(os.path.join(cfg.output_dir, f"{cfg.outcome}_metrics.json"), "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nDone. Outputs -> {cfg.output_dir}/")
    return {"results": results, "interactions": inter, "X": X, "y": y}


if __name__ == "__main__":
    config = Config()
    # Example overrides:
    # config.outcome = "any_metabolic"
    # config.score_method = "rank"
    run(config)
