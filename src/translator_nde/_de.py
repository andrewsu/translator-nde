"""Vendored from /home/asu/Science/DN-meta-analysis/scripts/lib.py.

Copied rather than imported so this repo is self-contained and its results stay
reproducible if that repo moves. Unmodified except for this header; see
tests/test_de_reproduces_dn.py, which checks it still reproduces the published
diabetic-nephropathy result.

Original docstring follows.
"""
"""Core library for the diabetic-nephropathy expression meta-analysis.

Everything is implemented on numpy/scipy so the pipeline is transparent and
does not depend on R/limma. Two statistical engines are provided:

  * limma-style empirical-Bayes moderated t-test  (per-study DE tables)
  * Hedges' g + DerSimonian-Laird random-effects   (cross-study pooling)
"""
import gzip, io, os, re, urllib.request
import numpy as np
import pandas as pd
from scipy import special, stats

# ----------------------------------------------------------------------------
# Download / parsing
# ----------------------------------------------------------------------------
def download(url, path, timeout=180):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    os.makedirs(os.path.dirname(path), exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "dn-meta/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r, open(path, "wb") as f:
        f.write(r.read())
    return path


def parse_series_matrix(path):
    """Return (expr[probe x GSM] float DataFrame, sample_meta dict[key->list])."""
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    sample_meta, table_lines, in_table = {}, [], False
    for line in text.splitlines():
        if line.startswith("!series_matrix_table_begin"):
            in_table = True
            continue
        if line.startswith("!series_matrix_table_end"):
            break
        if in_table:
            table_lines.append(line)
        elif line.startswith("!Sample_"):
            key = line.split("\t", 1)[0].strip("!")
            vals = [v.strip('"') for v in line.split("\t")[1:]]
            sample_meta.setdefault(key, [])
            # some keys (characteristics) repeat: store as list-of-lists
            if key in ("Sample_characteristics_ch1",):
                sample_meta[key].append(vals)
            else:
                sample_meta[key] = vals
    buf = io.StringIO("\n".join(table_lines))
    expr = pd.read_csv(buf, sep="\t", index_col=0, na_values=["", "null", "NA"])
    expr.index = expr.index.astype(str)
    expr = expr.apply(pd.to_numeric, errors="coerce")
    return expr, sample_meta


def maybe_log2(expr):
    """Log2-transform if values look linear (heuristic on 99th percentile)."""
    v = expr.values.astype(float)
    p99 = np.nanpercentile(v, 99)
    if p99 > 50:  # linear intensities
        shift = 0.0
        mn = np.nanmin(v)
        if mn <= 0:
            shift = 1.0 - mn
        return np.log2(expr + shift), True
    return expr, False


# ----------------------------------------------------------------------------
# Probe -> gene symbol mapping + collapse (max-mean probe per gene)
# ----------------------------------------------------------------------------
def parse_gene_assignment(s):
    """HTA-2.0 'gene_assignment' -> first gene symbol, or None."""
    if not isinstance(s, str) or "//" not in s:
        return None
    first = s.split("///")[0]
    parts = [p.strip() for p in first.split("//")]
    if len(parts) >= 2 and parts[1] and parts[1] != "---":
        return parts[1]
    return None


def build_probe2symbol(gpl_table, kind):
    """Return dict probe_id -> gene symbol for a platform annotation table."""
    t = gpl_table.copy()
    t["ID"] = t["ID"].astype(str)
    mapping = {}
    if kind == "gene_symbol_col":              # GPL571
        col = "Gene Symbol"
        for pid, sym in zip(t["ID"], t[col]):
            if isinstance(sym, str) and sym not in ("", "---"):
                mapping[pid] = sym.split("///")[0].strip()
    elif kind == "symbol_col":                 # GPL22945 (ENTREZG CDF)
        col = "Symbol"
        for pid, sym in zip(t["ID"], t[col]):
            if isinstance(sym, str) and sym not in ("", "---"):
                mapping[pid] = sym.strip()
    elif kind == "gene_assignment":            # GPL17586 (HTA 2.0)
        for pid, ga in zip(t["ID"], t["gene_assignment"]):
            sym = parse_gene_assignment(ga)
            if sym:
                mapping[pid] = sym
    else:
        raise ValueError(kind)
    return mapping


def collapse_by_symbol(expr):
    """Collapse a symbol-indexed matrix (possible duplicate symbols) to one row
    per gene, keeping the highest-mean row. Drops empty/invalid symbols."""
    idx = expr.index.astype(str)
    ok = [isinstance(s, str) and s not in ("", "---", "NA", "nan", "None") for s in idx]
    expr = expr.loc[ok].copy()
    expr.index = expr.index.astype(str)
    m = expr.mean(axis=1)
    expr = expr.assign(_m=m).sort_values("_m", ascending=False)
    expr = expr[~expr.index.duplicated(keep="first")].drop(columns="_m")
    expr.index.name = "gene"
    return expr.sort_index()


def counts_to_log2cpm(counts, min_cpm=1.0, min_frac=0.2):
    """Filter low-expressed genes then log2-CPM (limma-trend style)."""
    counts = counts.astype(float)
    lib_size = counts.sum(axis=0)
    cpm = counts.div(lib_size, axis=1) * 1e6
    keep = (cpm >= min_cpm).sum(axis=1) >= max(2, int(min_frac * counts.shape[1]))
    return np.log2(cpm.loc[keep] + 1.0)


# Immediate-early / stress-response genes: induced by tissue procurement, warm
# ischemia, and processing delay -> classic bulk-kidney artifact (dominated the
# ARCHS4 result). Flagged and excluded from the primary DEG set.
IEG_GENES = {
    "FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "EGR1", "EGR2",
    "EGR3", "EGR4", "NR4A1", "NR4A2", "NR4A3", "DUSP1", "DUSP2", "DUSP5",
    "DUSP6", "ATF3", "ARC", "IER2", "IER3", "IER5", "GADD45B", "GADD45G",
    "ZFP36", "BTG2", "CCN1", "CYR61", "CCN2", "CTGF", "RGS1", "RGS2", "PER1",
    "SOCS3", "KLF2", "KLF4", "KLF6", "MCL1", "HSPA1A", "HSPA1B", "DNAJB1",
    "NFKBIA", "PPP1R15A", "SGK1", "MYC", "JUND", "SIK1", "APOLD1", "MAFF",
}


def collapse_to_genes(expr, probe2sym):
    """Collapse probes to genes keeping, per gene, the probe with highest mean."""
    probes = expr.index.intersection(pd.Index(list(probe2sym.keys())))
    sub = expr.loc[probes]
    syms = pd.Series({p: probe2sym[p] for p in probes})
    means = sub.mean(axis=1)
    order = means.sort_values(ascending=False).index          # high-expressed first
    seen, keep_probe, keep_sym = set(), [], []
    for p in order:
        s = syms[p]
        if s not in seen:
            seen.add(s)
            keep_probe.append(p)
            keep_sym.append(s)
    genes = sub.loc[keep_probe].copy()
    genes.index = keep_sym
    genes.index.name = "gene"
    return genes.sort_index()


# ----------------------------------------------------------------------------
# limma-style empirical Bayes moderated t-test (two groups)
# ----------------------------------------------------------------------------
def _trigamma_inverse(x):
    """Solve trigamma(y)=x for y (Smyth 2004, fitFDist helper)."""
    x = np.asarray(x, float)
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        dif = tri * (1 - tri / x) / special.polygamma(2, y)
        y = y + dif
        if np.max(np.abs(dif / y)) < 1e-8:
            break
    return y


def fit_f_dist(s2, df):
    """Estimate prior (d0, s0^2) for gene variances (limma fitFDist)."""
    s2 = np.asarray(s2, float)
    ok = np.isfinite(s2) & (s2 > 0)
    s2 = s2[ok]
    z = np.log(s2)
    e = z - special.digamma(df / 2.0) + np.log(df / 2.0)
    ebar = e.mean()
    G = len(s2)
    evar = np.sum((e - ebar) ** 2) / (G - 1) - special.polygamma(1, df / 2.0)
    if evar > 0:
        d0 = 2 * _trigamma_inverse(evar)
        s0_2 = np.exp(ebar + special.digamma(d0 / 2.0) - np.log(d0 / 2.0))
    else:
        d0 = np.inf
        s0_2 = np.exp(ebar)
    return float(d0), float(s0_2)


def moderated_ttest(genes, is_dn):
    """limma-style two-group moderated t-test.  Positive logFC = up in DN.

    genes : DataFrame [gene x sample] (log2). is_dn : bool array over columns.
    Returns DataFrame with logFC, t, P.Value, adj.P.Val, AveExpr.
    """
    X = genes.values.astype(float)
    g1 = X[:, is_dn]                 # DN
    g2 = X[:, ~is_dn]                # control
    n1, n2 = g1.shape[1], g2.shape[1]
    m1, m2 = g1.mean(1), g2.mean(1)
    logfc = m1 - m2
    ss = ((g1 - m1[:, None]) ** 2).sum(1) + ((g2 - m2[:, None]) ** 2).sum(1)
    dg = n1 + n2 - 2
    s2 = ss / dg
    d0, s0_2 = fit_f_dist(s2, dg)
    if np.isinf(d0):
        s2_post, df_tot = s2, dg
    else:
        s2_post = (d0 * s0_2 + dg * s2) / (d0 + dg)
        df_tot = dg + d0
    unscaled = 1.0 / n1 + 1.0 / n2
    se = np.sqrt(s2_post * unscaled)
    t = logfc / se
    p = 2 * stats.t.sf(np.abs(t), df_tot)
    out = pd.DataFrame({
        "gene": genes.index, "AveExpr": X.mean(1), "logFC": logfc,
        "t": t, "P.Value": p, "n_dn": n1, "n_ctrl": n2,
    }).set_index("gene")
    out["adj.P.Val"] = bh_fdr(out["P.Value"].values)
    return out.sort_values("P.Value")


# ----------------------------------------------------------------------------
# Effect size (Hedges' g) per study
# ----------------------------------------------------------------------------
def hedges_g(genes, is_dn):
    X = genes.values.astype(float)
    g1, g2 = X[:, is_dn], X[:, ~is_dn]
    n1, n2 = g1.shape[1], g2.shape[1]
    m1, m2 = g1.mean(1), g2.mean(1)
    v1, v2 = g1.var(1, ddof=1), g2.var(1, ddof=1)
    sp = np.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
    sp = np.where(sp == 0, np.nan, sp)
    d = (m1 - m2) / sp
    J = 1 - 3.0 / (4 * (n1 + n2) - 9)
    g = J * d
    var_g = J ** 2 * ((n1 + n2) / (n1 * n2) + d ** 2 / (2 * (n1 + n2 - 2)))
    return pd.DataFrame({"gene": genes.index, "g": g, "var_g": var_g,
                         "logFC": m1 - m2, "n_dn": n1, "n_ctrl": n2}).set_index("gene")


# ----------------------------------------------------------------------------
# Meta-analysis
# ----------------------------------------------------------------------------
def bh_fdr(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full_like(p, np.nan)
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    ranked = np.clip(ranked, 0, 1)
    out = np.empty(n)
    out[order] = ranked
    q[ok] = out
    return q


def dersimonian_laird(g, v):
    """Random-effects pooling for one gene. g,v arrays over studies (finite)."""
    g, v = np.asarray(g, float), np.asarray(v, float)
    k = len(g)
    w = 1.0 / v
    g_fixed = np.sum(w * g) / np.sum(w)
    Q = np.sum(w * (g - g_fixed) ** 2)
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
    wr = 1.0 / (v + tau2)
    g_re = np.sum(wr * g) / np.sum(wr)
    se = np.sqrt(1.0 / np.sum(wr))
    z = g_re / se
    p = 2 * stats.norm.sf(abs(z))
    I2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
    return g_re, se, z, p, tau2, Q, I2, k


def meta_analyze(effect_tables, min_studies=2):
    """effect_tables: dict study_id -> DataFrame[gene]{g,var_g,logFC}.
    Returns per-gene random-effects meta result across studies."""
    studies = list(effect_tables)
    genes = set()
    for t in effect_tables.values():
        genes |= set(t.index)
    rows = []
    for gene in genes:
        gs, vs, lfc = [], [], []
        present = []
        for s in studies:
            t = effect_tables[s]
            if gene in t.index:
                gi, vi = t.at[gene, "g"], t.at[gene, "var_g"]
                if np.isfinite(gi) and np.isfinite(vi) and vi > 0:
                    gs.append(gi); vs.append(vi); lfc.append(t.at[gene, "logFC"])
                    present.append(s)
        if len(gs) < min_studies:
            continue
        g_re, se, z, p, tau2, Q, I2, k = dersimonian_laird(gs, vs)
        n_pos = int(np.sum(np.array(gs) > 0)); n_neg = int(np.sum(np.array(gs) < 0))
        rows.append({
            "gene": gene, "k": k, "pooled_g": g_re, "se": se, "z": z,
            "P.Value": p, "tau2": tau2, "Q": Q, "I2": I2,
            "mean_logFC": float(np.mean(lfc)),
            "n_up": n_pos, "n_down": n_neg,
            "dir_consistent": (n_pos == k) or (n_neg == k),
            "studies": ",".join(present),
        })
    res = pd.DataFrame(rows).set_index("gene")
    res["adj.P.Val"] = bh_fdr(res["P.Value"].values)
    return res.sort_values("P.Value")


def stouffer(pmat_signed, weights):
    """Secondary check: weighted Stouffer combination.
    pmat_signed: DataFrame[gene x study] of signed z-scores (NaN if absent)."""
    Z = pmat_signed.values
    W = np.array(weights, float)[None, :]
    mask = np.isfinite(Z)
    zsum = np.nansum(np.where(mask, Z * W, 0.0), axis=1)
    wnorm = np.sqrt(np.nansum(np.where(mask, W ** 2, 0.0), axis=1))
    zc = zsum / wnorm
    p = 2 * stats.norm.sf(np.abs(zc))
    return pd.DataFrame({"gene": pmat_signed.index, "stouffer_z": zc,
                         "stouffer_p": p}).set_index("gene")
