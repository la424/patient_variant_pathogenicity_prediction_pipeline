"""Output formatting: CSV and multi-sheet XLSX."""

import pandas as pd
from .config import PipelineConfig


def save_results(df: pd.DataFrame, config: PipelineConfig):
    """Save all output files."""
    out = config.output_settings
    results_dir = config.results_dir
    partner_labels = config.partner_labels

    # Column ordering
    id_cols = ['gene', 'position', 'ref_aa', 'alt_aa']
    grantham_cols = ['grantham_distance', 'grantham_class', 'substitution_severity', 'property_changes']
    summary_cols = [c for c in ['variant_summary', 'pipeline_agreement'] if c in df.columns]

    v6_score_cols = [c for c in [
        'v6_final_score', 'v6_tier', 'v6_contact_disruption', 'v6_total_contacts',
        'v6_max_inter_contacts', 'v6_best_inter_partner',
        'v6_best_burial', 'v6_best_burial_source',
        'v6_interface_partners_gated', 'v6_n_interface_partners_gated',
        'v6_score_evidence',
    ] if c in df.columns]

    mech_cols = [c for c in ['v6_mechanism', 'v6_mechanism_partner', 'v6_external_evidence_flag'] if c in df.columns]
    eval_cols = [c for c in ['structure_evaluable', 'ddg_evaluable', 'am_evaluable',
                             'franklin_evaluable', 'multimer_evaluable'] if c in df.columns]
    am_cols = [c for c in ['AlphaMissense', 'AlphaMissense_pathogenicity', 'AlphaMissense_raw'] if c in df.columns]
    franklin_cols = [c for c in ['franklin', 'franklin_raw'] if c in df.columns]
    gnomad_cols = [c for c in ['gnomad_af', 'gnomad_popmax_af', 'gnomad_homozygotes',
                               'gene_pli', 'gene_loeuf'] if c in df.columns]

    ddg_core_cols = [c for c in [
        'ddg_monomer', 'ddg_monomer_confident', 'ddg_category', 'ddg_category_raw', 'ddg_confidence',
        'ddg_multimer_max', 'ddg_multimer_min', 'ddg_multimer_mean',
        'ddg_category_multimer', 'ddg_category_multimer_raw',
        'n_complexes_tested', 'partners_tested',
        'ddg_low_confidence_flag',
        'ddg_summary_flags', 'ddg_summary_partners_affected', 'ddg_summary_partners_neutral',
    ] if c in df.columns]

    per_partner_cols = []
    for pl in partner_labels:
        for sfx in [f"ddg_binding_{pl}", f"ddg_fold_{pl}",
                    f"ddg_binding_category_{pl}", f"ddg_fold_category_{pl}",
                    f"ddg_{pl}_confident", f"ddg_interp_{pl}"]:
            if sfx in df.columns:
                per_partner_cols.append(sfx)

    conc_cols = [c for c in [
        'concordance_strict', 'concordance_relaxed', 'concordance_t3',
        'four_way_strict', 'four_way_strict_denom',
        'four_way_relaxed', 'four_way_relaxed_denom',
        'four_way_t3', 'four_way_t3_denom',
        'structure_vote_strict', 'structure_vote_t3',
        'ddg_vote_strict', 'ddg_vote_relaxed',
        'am_vote_strict', 'am_vote_relaxed',
        'franklin_vote_strict', 'franklin_vote_relaxed', 'max_abs_ddg',
    ] if c in df.columns]

    nbhd_cols = [c for c in df.columns if c.startswith('nbhd_')]
    mono_cols = [c for c in df.columns if c.startswith('monomer_')]
    multi_summary = [c for c in [
        'n_multimer_complexes', 'multimer_partners',
        'is_interface_any', 'interface_partners', 'n_interface_partners',
        'multimer_plddt_max', 'multimer_plddt_avg',
        'multimer_contacts_max', 'multimer_contacts_avg',
        'multimer_disruption_max', 'multimer_disruption_avg',
    ] if c in df.columns]
    multi_cols = [c for c in df.columns if c.startswith('multi_')]

    ordered = (
        id_cols + grantham_cols + summary_cols +
        v6_score_cols + mech_cols + eval_cols +
        am_cols + franklin_cols + gnomad_cols +
        ddg_core_cols + per_partner_cols +
        conc_cols + nbhd_cols + mono_cols + multi_summary + multi_cols
    )
    seen = set()
    final_ordered = []
    for c in ordered:
        if c in df.columns and c not in seen:
            final_ordered.append(c)
            seen.add(c)
    for c in df.columns:
        if c not in seen:
            final_ordered.append(c)
            seen.add(c)

    df_out = df[final_ordered].copy()

    # CSV
    csv_path = results_dir / out['csv_filename']
    df_out.to_csv(csv_path, index=False)
    print(f"✓ CSV: {csv_path} ({len(df_out)} rows, {len(df_out.columns)} columns)")

    # High priority
    hp = df_out[df_out['v6_tier'].isin(['Tier 1', 'Tier 2'])]
    hp_path = results_dir / out['high_priority_filename']
    hp.to_csv(hp_path, index=False)
    print(f"✓ High priority: {hp_path} ({len(hp)} variants)")

    # DDG expanded
    ddg_exp_cols = [c for c in [
        'gene', 'position', 'ref_aa', 'alt_aa', 'ddg_monomer', 'ddg_monomer_confident',
        'ddg_category', 'ddg_category_multimer',
        'ddg_multimer_max', 'ddg_multimer_min', 'ddg_multimer_mean',
        'n_complexes_tested', 'partners_tested',
        'ddg_summary_flags', 'ddg_summary_partners_affected',
        'ddg_low_confidence_flag',
    ] if c in df.columns]
    ddg_path = results_dir / out['ddg_expanded_filename']
    df[ddg_exp_cols].to_csv(ddg_path, index=False)
    print(f"✓ DDG expanded: {ddg_path}")

    # XLSX
    if out.get('xlsx_sheets', True):
        try:
            import openpyxl
            xlsx_path = results_dir / out['xlsx_filename']

            summary_sheet_cols = [c for c in (
                id_cols + ['grantham_distance'] +
                ['v6_final_score', 'v6_tier', 'v6_mechanism', 'v6_mechanism_partner'] +
                ['concordance_strict', 'concordance_relaxed', 'concordance_t3'] +
                ['AlphaMissense', 'franklin',
                 'ddg_monomer', 'ddg_monomer_confident', 'ddg_category',
                 'ddg_category_multimer', 'ddg_low_confidence_flag',
                 'ddg_summary_flags', 'ddg_summary_partners_affected'] +
                ['structure_evaluable', 'ddg_evaluable', 'multimer_evaluable'] +
                gnomad_cols +
                ['nbhd_tier', 'nbhd_mechanism', 'nbhd_concordance_strict'] +
                ['pipeline_agreement', 'variant_summary', 'v6_external_evidence_flag']
            ) if c in df.columns]

            with pd.ExcelWriter(xlsx_path, engine='openpyxl') as writer:
                df[summary_sheet_cols].to_excel(writer, sheet_name='Summary', index=False)
                df[id_cols + [c for c in mono_cols if c in df.columns]].to_excel(
                    writer, sheet_name='Monomer Detail', index=False)
                df[id_cols + [c for c in (multi_summary + multi_cols) if c in df.columns]].to_excel(
                    writer, sheet_name='Multimer Detail', index=False)
                df[id_cols + [c for c in (ddg_core_cols + per_partner_cols) if c in df.columns]].to_excel(
                    writer, sheet_name='DDG Detail', index=False)
                df[id_cols + [c for c in (eval_cols + conc_cols + nbhd_cols) if c in df.columns]].to_excel(
                    writer, sheet_name='Concordance Detail', index=False)
            print(f"✓ XLSX: {xlsx_path} (5 sheets)")
        except ImportError:
            print("⚠ openpyxl not installed — XLSX skipped")

    # Summary
    print(f"\n{'=' * 60}")
    print(f"Pipeline complete — {len(df_out)} variants, {len(df_out.columns)} columns")
    print(f"{'=' * 60}")
    print(f"  Tiers: {df_out['v6_tier'].value_counts().to_dict()}")
    print(f"  Mechanisms: {df_out['v6_mechanism'].nunique()} categories")
    print(f"  Strict 4/4: {int((df_out.get('four_way_strict', pd.Series()) == 4).sum())}")
    print(f"{'=' * 60}")

    return df_out
