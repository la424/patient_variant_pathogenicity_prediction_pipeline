"""Main pipeline orchestrator."""

import warnings
import pandas as pd

from .config import PipelineConfig
from .extraction import extract_monomer_features, extract_multimer_features
from .ddg import (
    load_monomer_ddg, load_multimer_ddg, load_per_partner_ddg,
    fix_multimer_min_nan, normalize_partners_tested,
    compute_monomer_confidence, compute_per_partner_confidence,
    compute_ddg_summary, compute_multimer_category,
    compute_low_confidence_flags,
    run_missing_monomer_ddg, run_missing_multimer_ddg,
)
from .scoring import (
    compute_score, assign_tier, compute_evaluability,
    classify_mechanism, compute_concordance, normalize_annotations,
)
from .output import save_results


def run_pipeline(config_path: str) -> pd.DataFrame:
    """Run the complete variant pathogenicity pipeline.

    Args:
        config_path: Path to the YAML configuration file.

    Returns:
        Complete DataFrame with all pipeline columns.
    """
    warnings.filterwarnings('ignore')
    config = PipelineConfig(config_path)

    print("=" * 60)
    print("Variant Pathogenicity Pipeline v1.0")
    print("=" * 60)

    # ── Stage 1: Load variants ────────────────────────────────────────────
    print("\n[1/8] Loading variants...")
    df = pd.read_csv(config.variants_file)

    # Apply column mapping
    cm = config.column_map
    rename = {}
    for internal, external in cm.items():
        if external in df.columns and external != internal:
            rename[external] = internal
    if rename:
        df = df.rename(columns=rename)

    df.columns = [c.lower().strip() for c in df.columns]
    df['position'] = df['position'].astype(int)

    # Deduplicate
    n_before = len(df)
    df = df.drop_duplicates(subset=['gene', 'position', 'ref_aa', 'alt_aa']).reset_index(drop=True)
    if len(df) < n_before:
        print(f"  ⚠ Removed {n_before - len(df)} duplicate rows")

    # Preserve annotation columns
    ann_cols = [c for c in df.columns if c in [
        'alphamissense', 'franklin', 'alphamissense_pathogenicity',
        'gnomad_af', 'gnomad_popmax_af', 'gnomad_homozygotes',
        'gene_pli', 'gene_loeuf',
    ]]
    annotation_df = df[['gene', 'position', 'ref_aa', 'alt_aa'] + ann_cols].copy()

    print(f"  ✓ {len(df)} variants across {df['gene'].nunique()} genes")

    # ── Stage 2: Monomer extraction ───────────────────────────────────────
    print("\n[2/8] Extracting monomer features...")
    df = extract_monomer_features(df, config)

    # ── Stage 3: Multimer extraction ──────────────────────────────────────
    print("\n[3/8] Extracting multimer features...")
    df = extract_multimer_features(df, config)

    # ── Stage 4: Load DDG data ────────────────────────────────────────────
    print("\n[4/8] Loading DDG data...")
    df = load_monomer_ddg(df, config)
    df = load_multimer_ddg(df, config)

    # Ensure required columns exist
    for c in ['ddg_monomer', 'ddg_multimer_max', 'ddg_multimer_min',
              'ddg_multimer_mean', 'n_complexes_tested', 'partners_tested']:
        if c not in df.columns:
            df[c] = None

    df = fix_multimer_min_nan(df)
    df = normalize_partners_tested(df, config)

    # Monomer confidence
    df = compute_monomer_confidence(df, config)

    # Per-partner DDG
    df = load_per_partner_ddg(df, config)
    df = compute_per_partner_confidence(df, config)

    # ── Stage 5: Run missing FoldX (optional) ─────────────────────────────
    print("\n[5/8] Computing missing DDG values...")
    df = run_missing_monomer_ddg(df, config)
    df = run_missing_multimer_ddg(df, config)

    # Recompute confidence after new DDG values
    df = compute_monomer_confidence(df, config)
    df = compute_per_partner_confidence(df, config)

    # ── Stage 6: DDG summary and categories ───────────────────────────────
    print("\n[6/8] Computing DDG summaries...")
    df = compute_ddg_summary(df, config)
    df = compute_multimer_category(df, config)
    df = compute_low_confidence_flags(df, config)

    # ── Stage 7: Scoring, mechanism, concordance ──────────────────────────
    print("\n[7/8] Scoring and classification...")

    # Scoring
    score_results = df.apply(lambda r: compute_score(r, config), axis=1)
    for c in score_results.columns:
        df[c] = score_results[c]
    df['v6_tier'] = df['v6_final_score'].apply(lambda s: assign_tier(s, config))

    print(f"  Tiers: {df['v6_tier'].value_counts().to_dict()}")

    # Evaluability
    df = compute_evaluability(df, config)

    # Mechanism
    mech_results = df.apply(lambda r: classify_mechanism(r, config), axis=1)
    df['v6_mechanism'] = [r[0] for r in mech_results]
    df['v6_mechanism_partner'] = [r[1] for r in mech_results]
    df['v6_external_evidence_flag'] = [r[2] for r in mech_results]

    print(f"  Mechanisms:\n{df['v6_mechanism'].value_counts().to_string()}")

    # Merge annotations and normalize
    df = _merge_annotations(df, annotation_df)
    df = normalize_annotations(df, config)

    # Concordance
    conc_results = df.apply(lambda r: compute_concordance(r, config), axis=1)
    for c in conc_results.columns:
        df[c] = conc_results[c]

    strict_44 = int((df.get('four_way_strict', pd.Series()) == 4).sum())
    print(f"  Strict 4/4: {strict_44}")

    # Variant summary string
    df['variant_summary'] = df.apply(_variant_summary, axis=1)

    # ── Stage 8: Save ─────────────────────────────────────────────────────
    print("\n[8/8] Saving results...")
    df_out = save_results(df, config)

    return df_out


def _merge_annotations(df, annotation_df):
    """Merge annotation columns back from original data."""
    # Drop any existing annotation cols to prevent duplication
    drop_cols = [c for c in ['AlphaMissense', 'AlphaMissense_pathogenicity',
                             'AlphaMissense_raw', 'franklin', 'franklin_raw',
                             'alphamissense', 'alphamissense_pathogenicity']
                 if c in df.columns]
    if drop_cols:
        df = df.drop(columns=drop_cols)

    df['gene_lower'] = df['gene'].astype(str).str.lower()
    ann = annotation_df.copy()
    ann.columns = [c.lower() for c in ann.columns]
    ann['gene_lower'] = ann['gene'].astype(str).str.lower()

    for col in ['alphamissense', 'alphamissense_pathogenicity', 'franklin']:
        if col in ann.columns:
            m = ann[['gene_lower', 'position', 'ref_aa', 'alt_aa', col]].drop_duplicates()
            df = df.merge(m, on=['gene_lower', 'position', 'ref_aa', 'alt_aa'], how='left')
            rn = {'alphamissense': 'AlphaMissense',
                  'alphamissense_pathogenicity': 'AlphaMissense_pathogenicity'}
            if col in rn:
                df = df.rename(columns={col: rn[col]})

    df = df.drop(columns=['gene_lower'], errors='ignore')

    for col in ['AlphaMissense', 'AlphaMissense_pathogenicity', 'franklin']:
        if col not in df.columns:
            df[col] = None
    return df


def _variant_summary(row):
    """Generate human-readable variant summary string."""
    parts = []
    tier = row.get('v6_tier', '')
    mech = row.get('v6_mechanism', '')
    conc = row.get('concordance_strict', '')
    am = str(row.get('AlphaMissense', ''))
    fr = str(row.get('franklin', ''))

    parts.append(f"{tier}")
    if mech and mech != 'Structure unevaluable':
        parts.append(mech)
    if conc:
        parts.append(f"Concordance:{conc}")
    if am and am != 'nan':
        parts.append(f"AM:{am}")
    if fr and fr != 'nan' and fr != 'No data':
        parts.append(f"Franklin:{fr}")

    flags = str(row.get('ddg_summary_flags', ''))
    if flags and 'fold_unknown' not in flags:
        parts.append(flags.split(';')[0])

    return ' | '.join(parts)
