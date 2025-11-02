# ------------------------------------------------------------
# AUTHORSHIP
# ------------------------------------------------------------
# Author: Hashim Al-Hashimi Lab
# Created: October 2025
# ------------------------------------------------------------

# ------------------------------------------------------------
# REQUIREMENTS NOTICE
# ------------------------------------------------------------
# This script requires Python 3.8 or higher and the following
# package versions (or newer) to ensure full compatibility:
#
#   - numpy >= 1.20
#   - pandas >= 1.1
#   - scipy >= 1.8        
#   - matplotlib >= 3.3
# ------------------------------------------------------------

import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon
from scipy.interpolate import UnivariateSpline

# -------------------------------------------
# Visualization settings
# -------------------------------------------
CANONICAL_ORDER = [
    "ATA", "ATC", "ATG", "ATT",
    "CTA", "CTC", "CTG", "CTT",
    "GTA", "GTC", "GTG", "GTT",
    "TTA", "TTC", "TTG", "TTT"
]

# -------------------------------------------
# File I/O
# -------------------------------------------
def load_data(conf_path, mut_sig_path):
    params_conf = pd.read_csv(
        conf_path, sep=",", skiprows=1, names=["Trinucleotide", "pB"]
    )
    params_conf["pB"] = pd.to_numeric(params_conf["pB"], errors="coerce")
    params_conf = params_conf[["Trinucleotide", "pB"]]
    mut_sig_profiles = pd.read_csv(mut_sig_path, sep="\t")
    return params_conf, mut_sig_profiles

# -------------------------------------------
# Convert purine-centered → pyrimidine-centered
# -------------------------------------------
def conv_seq_A_to_T(seq_con_A):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    # seq_con_A is an iterable of strings like 'ATA' etc.
    return [complement.get(s[2], "?") + "T" + complement.get(s[0], "?") for s in seq_con_A]

# -------------------------------------------
# Process COSMIC signatures into substitution groups
# -------------------------------------------
def process_signatures(mut_sig_profiles):
    subs_order = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    parsed = mut_sig_profiles["Type"].str.extract(r'([ACGT])\[(.+)\]([ACGT])')
    mut_sig_profiles = mut_sig_profiles.copy()
    mut_sig_profiles["sub_type"] = parsed[1]
    mut_sig_profiles["seqs_mut"] = parsed[0] + "T" + parsed[2]
    mut_sig_profiles.sort_values(
        by="sub_type",
        key=lambda x: x.map({s: i for i, s in enumerate(subs_order)}),
        inplace=True
    )
    # Return dictionary mapping substitution string -> dataframe (with seqs_mut & signature columns)
    return {sub: mut_sig_profiles.loc[mut_sig_profiles["sub_type"] == sub].copy()
            for sub in subs_order}

# -------------------------------------------
# Utils
# -------------------------------------------
def normalize_keep_shape(values, mask):
    """
    values: 1D np array-like
    mask: boolean array (same length) indicating which entries should be kept/normalized
    Returns array same length: non-selected entries set to 0; selected entries normalized to sum to 1.
    """
    arr = np.array(values, dtype=float)
    arr[~mask] = 0.0
    if mask.sum() > 0:
        s = arr[mask].sum()
        if s > 0:
            arr[mask] = arr[mask] / s
        else:
            arr[mask] = 0.0
    return arr

def interpret_sub_filter(filter_value):
    full = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    if filter_value is None:
        return full
    if isinstance(filter_value, str):
        s = filter_value.strip().upper()
        if s in full:
            return [s]
        raise ValueError(f"Unrecognized filter string: {filter_value}")
    if isinstance(filter_value, (list, tuple)):
        result = []
        for x in filter_value:
            sx = str(x).upper()
            if sx in full:
                result.append(sx)
            else:
                raise ValueError(f"Unrecognized filter string in list: {x}")
        return result
    raise ValueError("filter_value must be None, str, or list/tuple")

# -------------------------------------------
# Build SBS table for a (possibly combined) substitution set
# -------------------------------------------
def _combined_df_for_sub(grouped, sub):
    """
    Return (df_sub, sig_cols) where df_sub has columns ['seqs_mut'] + sig_cols.
    """
    if sub not in grouped:
        return None, []
    df = grouped[sub]
    non_sig_cols = {"Type", "sub_type", "seqs_mut"}
    sig_cols = [c for c in df.columns if c not in non_sig_cols]
    return df[["seqs_mut"] + sig_cols].copy(), sig_cols

# -------------------------------------------
# Compute JSDs (observed)
# -------------------------------------------
def compute_js_divergence(grouped, params_conf, sub_filter=None, pyrimidine_centered=False):
    """
    Compute observed JSDs for all requested substitution sets.
    Returns dict mapping substitution -> {"sig_names": [...], "js": np.array([...])}
    """
    subs_to_run = interpret_sub_filter(sub_filter)
    js_div = {}

    params_conf = params_conf.copy()
    if not pyrimidine_centered:
        params_conf["Trinucleotide"] = conv_seq_A_to_T(params_conf["Trinucleotide"])

    for sub in subs_to_run:
        df_sub, sig_cols = _combined_df_for_sub(grouped, sub)
        if df_sub is None:
            js_div[sub] = {"sig_names": [], "js": np.array([])}
            continue

        # Merge to keep only contexts present in both
        merged = pd.merge(
            df_sub, params_conf,
            left_on="seqs_mut", right_on="Trinucleotide", how="inner"
        )

        js_values = []
        for sig in sig_cols:
            P_raw = pd.to_numeric(merged[sig], errors="coerce").values
            Q_raw = pd.to_numeric(merged["pB"], errors="coerce").values
            mask = (~np.isnan(P_raw)) & (~np.isnan(Q_raw))
            if mask.sum() == 0:
                js_values.append(np.nan)
                continue
            # Build masked normalized vectors for comparison (length = mask.sum())
            P_mask = P_raw[mask].astype(float)
            Q_mask = Q_raw[mask].astype(float)
            if P_mask.sum() <= 0 or np.all(np.isnan(Q_mask)):
                js_values.append(np.nan)
                continue
            P_norm = P_mask / P_mask.sum()
            Q_norm = Q_mask / Q_mask.sum()
            js_values.append((jensenshannon(P_norm, Q_norm, base=2))**2)

        js_div[sub] = {"sig_names": sig_cols, "js": np.array(js_values)}

    return js_div

# -------------------------------------------
# Helper: Storey q-values (λ-grid + spline smoothing)
# -------------------------------------------
def storey_qvalues(pvals, lambda_grid=None):
    """
    Compute Storey q-values using spline smoothing over lambda_grid.
      - lambda_grid defaults to 0,0.05,...,0.95
      - compute pi0_hat(lambda) = mean(p > lambda) / (1 - lambda)
      - fit cubic spline to pi0_hat vs lambda and evaluate at max(lambda)
      - q = pi0 * m * p / rank
      - enforce monotonicity
    """
    pvals = np.asarray(pvals, dtype=float)
    if lambda_grid is None:
        lambda_grid = np.arange(0.0, 0.95 + 1e-12, 0.05)  # 0.00..0.95 inclusive

    m = float(len(pvals))
    pi0_vals = []
    for l in lambda_grid:
        pi0_l = np.mean(pvals > l) / (1.0 - l) if (1.0 - l) > 0 else 1.0
        pi0_vals.append(pi0_l)
    pi0_vals = np.clip(np.array(pi0_vals, dtype=float), 0.0, 1.0)

    # Smooth using a cubic spline
    try:
        spline = UnivariateSpline(lambda_grid, pi0_vals, k=3, s=None)
        pi0 = float(spline(lambda_grid[-1]))
    except Exception:
        pi0 = float(pi0_vals[-1])

    pi0 = max(0.0, min(1.0, pi0))
    ranks = pd.Series(pvals).rank(method="average").values
    q_vals = (pi0 * m * pvals) / ranks
    order_desc = np.argsort(-pvals)
    q_desc = q_vals[order_desc]
    q_desc_mon = np.minimum.accumulate(q_desc)
    q_vals[order_desc] = q_desc_mon
    q_vals = np.clip(q_vals, 0.0, 1.0)
    return q_vals

# -------------------------------------------
# Helper: BH q-values with average ranks + monotonicity
# -------------------------------------------
def bh_qvalues(pvals):
    """
    Compute standard BH q-values:
      - average ranks for ties
      - q = p * m / rank
      - enforce monotonicity
    """
    pvals = np.asarray(pvals, dtype=float)
    m = float(len(pvals))
    ranks = pd.Series(pvals).rank(method="average").values
    q_vals = (pvals * m) / ranks
    order_desc = np.argsort(-pvals)
    q_desc = q_vals[order_desc]
    q_desc_mon = np.minimum.accumulate(q_desc)
    q_vals[order_desc] = q_desc_mon
    q_vals = np.clip(q_vals, 0.0, 1.0)
    return q_vals

# -------------------------------------------
# Permutation test + FDR (BH or Storey)
# -------------------------------------------
def permutation_test_and_fdr(js_div, grouped, params_conf, n_permutations=10000, 
                                         fdr_alpha=0.05, pyrimidine_centered=False, seed=1, 
                                         fdr_method="storey"):
    
    rng = np.random.default_rng(seed)
    params_conf_local = params_conf.copy()

    if not pyrimidine_centered:
        params_conf_local["Trinucleotide"] = conv_seq_A_to_T(params_conf_local["Trinucleotide"])

    params_conf_local = params_conf_local.set_index("Trinucleotide").reindex(CANONICAL_ORDER).reset_index()
    params_conf_local["Trinucleotide"] = CANONICAL_ORDER

    subs_order = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]

    mut_cols, mut_col_names = [], []
    for sub in subs_order:
        if sub not in grouped:
            continue
        df_sub = grouped[sub].set_index("seqs_mut")
        sig_cols = [c for c in df_sub.columns if c not in {"Type", "sub_type", "seqs_mut"}]
        df_reindexed = df_sub.reindex(CANONICAL_ORDER)
        for sig in sig_cols:
            mut_cols.append(df_reindexed[sig].to_numpy())
            mut_col_names.append(f"{sub}::{sig}")

    if not mut_cols:
        raise RuntimeError("No signature columns found in grouped; cannot build mutational matrix.")

    mut_matrix = np.column_stack(mut_cols)
    col_sums = np.nansum(mut_matrix, axis=0, keepdims=True)
    mut_matrix = np.divide(mut_matrix, col_sums, where=col_sums != 0)

    row_finite_values = [row[np.isfinite(row)] for row in mut_matrix]
    global_null_js = np.empty(n_permutations, dtype=float)

    Q_raw = params_conf_local["pB"].to_numpy(dtype=float)
    Q_raw = Q_raw / np.nansum(Q_raw)
    Q_mask = np.isfinite(Q_raw)
    Q_base = Q_raw[Q_mask] / np.sum(Q_raw[Q_mask])

    n_contexts = len(row_finite_values)
    mut_sig_rand_mat = np.zeros((n_contexts, n_permutations), dtype=float)

    for j, finite_vals in enumerate(row_finite_values):
        if finite_vals.size > 0:
            mut_sig_rand_mat[j, :] = rng.choice(finite_vals, size=n_permutations)
        else:
            mut_sig_rand_mat[j, :] = 0.0

    mut_sig_rand_mat /= np.sum(mut_sig_rand_mat, axis=0, keepdims=True)
    global_null_js = np.empty(n_permutations, dtype=float)

    for x in range(n_permutations):
        global_null_js[x] = jensenshannon(mut_sig_rand_mat[Q_mask, x], Q_base, base=2)**2

    global_null_js = global_null_js[np.isfinite(global_null_js)]
    global_null_js.sort()

    print(f"Global null JSDs: mean={global_null_js.mean():.5f}, "
          f"std={global_null_js.std():.5f}, "
          f"min={global_null_js.min():.5f}, max={global_null_js.max():.5f}")
    
    results = []

    for sub, info in js_div.items():
        sig_names, js_vals = info["sig_names"], info["js"]
        if not sig_names:
            continue

        df_sub, _ = _combined_df_for_sub(grouped, sub)
        if df_sub is None:
            continue

        for name, obs_js in zip(sig_names, js_vals):
            if np.isnan(obs_js):
                continue

            Count_Leq = int(np.sum(global_null_js <= obs_js))
            p_val = Count_Leq / global_null_js.size

            results.append([
                name, sub, float(obs_js), p_val,
                int(global_null_js.size), Count_Leq
            ])

    df = pd.DataFrame(results, columns=[
        "Signature", "Substitution_Type", "JSD", "p_value",
        "Null_N", "Count_Leq"
    ])

    df = df.dropna(subset=["p_value"]).reset_index(drop=True)
    
    if not df.empty:
        pvals = df["p_value"].to_numpy()
        if fdr_method.lower() == "bh":
            q_vals = bh_qvalues(pvals)
        elif fdr_method.lower() == "storey":
            q_vals = storey_qvalues(pvals)
        else:
            raise ValueError("fdr_method must be 'bh' or 'storey'")

        df["q_value"] = q_vals
        df["Pass_FDR"] = df["q_value"] < fdr_alpha
    else:
        df["q_value"] = []
        df["Pass_FDR"] = []

    return df

# -------------------------------------------
# Export observed JSD summary (CSV)
# -------------------------------------------
def export_js_summary(results_df, output_csv="JSD_summary_unfiltered.csv"):
    sub_order = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    if "Substitution_Type" in results_df.columns:
        results_df["Substitution_Type"] = pd.Categorical(
            results_df["Substitution_Type"], categories=sub_order, ordered=True
        )
        results_df = results_df.sort_values(by=["Substitution_Type", "JSD"]).reset_index(drop=True)
        results_df["Rank"] = results_df.groupby("Substitution_Type", observed=False)["JSD"].rank(method="dense")
    
    for col in ["p_value", "q_value"]:
        if col in results_df.columns:
            results_df[col] = results_df[col].astype(float)

    results_df.to_csv(output_csv, index=False, float_format="%.8f")
    print(f"✅ JSD summaries CSV saved to: {os.path.abspath(output_csv)}")
    print("Top lines of summary CSV:")
    print(results_df.head(10))
    return results_df

def export_filtered_js_summary(results_df, js_threshold, fdr_alpha, output_csv="JSD_summary_filtered.csv"):
    """
    Save subset of results that pass both JSD threshold and FDR q-value cutoff,
    grouped by substitution type.
    """

    if js_threshold is None:
        js_mask = np.ones(len(results_df), dtype=bool)
    else:
        js_mask = results_df["JSD"] < js_threshold

    filtered = results_df[(js_mask) & (results_df["q_value"] < fdr_alpha)].copy()
    if filtered.empty:
        print("⚠️ No results passed both JSD and FDR thresholds. No filtered CSV saved.")
        return filtered

    for col in ["p_value", "q_value"]:
        if col in filtered.columns:
            filtered[col] = filtered[col].astype(float)

    if "Substitution_Type" in filtered.columns:
        filtered = filtered.sort_values(
            by=["Substitution_Type", "q_value", "JSD"],
            ascending=[True, True, True]
        )

    filtered.to_csv(output_csv, index=False, float_format="%.8f")

    print(f"✅ Filtered JSD summary CSV saved to: {os.path.abspath(output_csv)}")
    print("Top lines of filtered summary CSV (grouped by Substitution_Type):")
    print(filtered.head(10))
    return filtered

# -------------------------------------------
# Main
# -------------------------------------------
def main():
    # ------------------------------------------------------------
    # INPUT FILE EXPECTATIONS
    # ------------------------------------------------------------
    # This script expects two input files:
    #
    # 1. Background probability file (pB file)
    #    - Format: CSV or tab-delimited text file (.csv, .tsv, .txt)
    #    - Must include:
    #        "Trinucleotide"  -> sequence context identifiers
    #        "pB"              -> background probabilities
    #    - The delimiter in conf_path should be specified in 
    #      the load_data function (sep=";" or sep="," or sep="/t")
    #    - If a trinucleotide context is missing a pB value, set pB to
    #      be "NaN" (e.g., "AAA;NaN")
    #    - Example: input_background.csv or input_background.tsv
    #
    # 2. COSMIC mutational signature database
    #    - Format: CSV or text file (.csv, .txt)
    #    - Must include:
    #        "Type"            -> substitution context in the form A[C>T]G
    #    - All other columns should correspond to signature probability profiles
    #      (e.g., SBS1, SBS2, SBS3, ...).
    #    - https://cancer.sanger.ac.uk/signatures/downloads/
    #    - https://doi.org/10.1038/s41586-020-1943-3
    #
    # Unsupported formats:
    #    - Excel files (.xlsx, .xls) and JSON or Parquet files are not read automatically.
    #    - Missing "Type" or "Trinucleotide" columns will cause KeyErrors.
    #
    # Notes:
    #    - Files containing NaN values are handled automatically.
    #    - Use UTF-8 encoding when saving files to avoid UnicodeDecodeError.
    # ------------------------------------------------------------
    
    parent_path = "/Users/orsula1/OrsDocs/OrsDocs_since_Jan2025/2025-Sequence-dependence/2025-Data-Deposition/Mutational-Signatures-JSD/JSD-python/"
    
    conf_input_path = "Input-files/conf-sig-input-files-for-python-JSDs/"
    mut_sig_input_path = "Input-files/mut-sig-input-files-for-python-JSDs/"
    
    conf_filename = "GT_Anion.csv"
    mut_sig_filename = "GRCh37_SBS_profiles.txt"
    
    conf_path = os.path.join(parent_path,conf_input_path,conf_filename)
    mut_sig_path = os.path.join(parent_path,mut_sig_input_path,mut_sig_filename)
    
    # OUTPUT Folder
    jsd_output_path = "Output-JSDs/GT-anion-JSDs/"
    jsd_output_filename_filtered = "JSD_summary_filtered.csv"
    jsd_output_filename_unfiltered = "JSD_summary_unfiltered.csv"
    
    jsd_f_path = os.path.join(parent_path,jsd_output_path,jsd_output_filename_filtered)
    jsd_uf_path = os.path.join(parent_path,jsd_output_path,jsd_output_filename_unfiltered)

    # Substitution filter
    # Examples:
    # sub_filter = "T>A"
    # sub_filter = None    # all six separately
    sub_filter = None

    js_plot_threshold = 0.090   # Threshold for the JSD test. None to disable
    fdr_alpha = 0.05            # FDR significance level
    fdr_method = "bh"           # Use "storey" for storey FDR test or "bh" for standard FDR test
    n_permutations = 1000000    # Number of permutations for null distribution.

    # Flag to indicate whether input sequences are already pyrimidine-centered (T or C)
    pyrimidine_centered = True  # Set to True if input sequences are already T/C centered

    # ----------------------
    # Load data, process signatures, compute JSD, run permutation/FDR, export summary, and plot
    # ----------------------
    params_conf, mut_sig_profiles = load_data(conf_path, mut_sig_path)
    grouped = process_signatures(mut_sig_profiles)

    # Observed JSDs
    js_div = compute_js_divergence(grouped, params_conf, sub_filter=sub_filter,
                                   pyrimidine_centered=pyrimidine_centered)

    # Permutation test + FDR (per-signature nulls via permuting pB)
    results_df = permutation_test_and_fdr(js_div, grouped, params_conf,
                                          n_permutations=n_permutations,
                                          fdr_alpha=fdr_alpha,
                                          pyrimidine_centered=pyrimidine_centered,
                                          seed=1,
                                          fdr_method=fdr_method)

    # Export CSV with JSD, p-value, q-value and rank
    export_js_summary(results_df, jsd_uf_path)
    export_filtered_js_summary(results_df, js_plot_threshold, fdr_alpha, jsd_f_path)

if __name__ == "__main__":
    main()
