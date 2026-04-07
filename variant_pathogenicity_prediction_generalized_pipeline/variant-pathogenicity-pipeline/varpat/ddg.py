"""DDG loading, three-axis framework, per-partner processing, and summary computation."""

import pandas as pd
import numpy as np
from pathlib import Path

from .config import PipelineConfig
from .constants import classify_ddg, sf
from .foldx import run_buildmodel, run_multimer_ddg


def load_monomer_ddg(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    """Load pre-computed monomer DDG values."""
    ddg_file = config.monomer_ddg_file
    if ddg_file is None or not ddg_file.exists():
        if 'ddg_monomer' not in df.columns:
            df['ddg_monomer'] = None
        print(f"  Monomer DDG file not found (existing: {df['ddg_monomer'].notna().sum()}/{len(df)})")
        return df

    ddg = pd.read_csv(ddg_file)
    ddg.columns = [c.lower() for c in ddg.columns]
    dc = next((c for c in ['ddg', 'ddg_monomer', 'total_ddg'] if c in ddg.columns), None)
    if dc:
        ddg = ddg.rename(columns={dc: 'ddg_monomer'})
        if 'ddg_monomer' in df.columns:
            df = df.drop(columns=['ddg_monomer'])
        df = df.merge(ddg[['gene', 'position', 'ref_aa', 'alt_aa', 'ddg_monomer']],
                      on=['gene', 'position', 'ref_aa', 'alt_aa'], how='left')
        print(f"✓ Monomer DDG: {df['ddg_monomer'].notna().sum()}/{len(df)}")
    return df


def load_multimer_ddg(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    """Load pre-computed aggregate multimer DDG values."""
    ddg_file = config.multimer_ddg_file
    if ddg_file is None or not ddg_file.exists():
        for c in ['ddg_multimer_max', 'ddg_multimer_min', 'ddg_multimer_mean',
                  'n_complexes_tested', 'partners_tested']:
            if c not in df.columns:
                df[c] = None
        return df

    ddgm = pd.read_csv(ddg_file)
    ddgm.columns = [c.lower() for c in ddgm.columns]
    dc = next((c for c in ['ddg', 'ddg_multimer', 'total_ddg'] if c in ddgm.columns), None)
    if dc:
        grp = ddgm.groupby(['gene', 'position', 'ref_aa', 'alt_aa'])
        agg = grp.agg(
            ddg_multimer_max=(dc, 'max'), ddg_multimer_min=(dc, 'min'),
            ddg_multimer_mean=(dc, 'mean'), n_complexes_tested=(dc, 'count')
        ).reset_index()
        if 'partner' in ddgm.columns:
            pt = grp['partner'].apply(lambda x: ';'.join(x.astype(str))).reset_index()
            pt.columns = ['gene', 'position', 'ref_aa', 'alt_aa', 'partners_tested']
            agg = agg.merge(pt, on=['gene', 'position', 'ref_aa', 'alt_aa'], how='left')
        for c in agg.columns:
            if c in df.columns and c not in ['gene', 'position', 'ref_aa', 'alt_aa']:
                df = df.drop(columns=[c])
        df = df.merge(agg, on=['gene', 'position', 'ref_aa', 'alt_aa'], how='left')
        print(f"✓ Multimer DDG: {df['ddg_multimer_max'].notna().sum()}/{len(df)}")
    return df


def load_per_partner_ddg(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    """Load per-partner DDG (binding + fold) from CSV."""
    ddg_file = config.per_partner_ddg_file
    if ddg_file is None or not ddg_file.exists():
        print("  No per-partner DDG file found")
        return df

    ppd = pd.read_csv(ddg_file)
    ppd.columns = [c.lower() for c in ppd.columns]

    # Normalize partner names
    if 'partner' in ppd.columns:
        ppd['partner'] = ppd.apply(
            lambda r: config.normalize_partner(r['partner'], r['gene']), axis=1)

    n_loaded = 0
    for _, row in ppd.iterrows():
        gene = str(row['gene']).lower()
        pos = int(row['position'])
        partner = str(row.get('partner', '')).lower()
        if not partner:
            continue

        mask = (df['gene'].str.lower() == gene) & (df['position'] == pos)
        if mask.sum() == 0:
            continue

        idx = df.loc[mask].index[0]
        for src_col, tgt_sfx in [('ddg_binding', 'ddg_binding'),
                                  ('ddg_fold', 'ddg_fold'),
                                  ('ddg', 'ddg_binding')]:  # fallback
            if src_col in row.index and pd.notna(row[src_col]):
                df.at[idx, f"{tgt_sfx}_{partner}"] = float(row[src_col])
                n_loaded += 1

    print(f"✓ Per-partner DDG: {n_loaded} values loaded")
    return df


def fix_multimer_min_nan(df: pd.DataFrame) -> pd.DataFrame:
    """Reconstruct missing ddg_multimer_min from max and mean."""
    n_fixed = 0
    for idx in df.index:
        if pd.notna(df.at[idx, 'ddg_multimer_max']) and pd.isna(df.at[idx, 'ddg_multimer_min']):
            n = df.at[idx, 'n_complexes_tested']
            mx = float(df.at[idx, 'ddg_multimer_max'])
            mn = df.at[idx, 'ddg_multimer_mean']
            if pd.notna(n) and pd.notna(mn):
                n, mn = int(n), float(mn)
                if n == 1:
                    df.at[idx, 'ddg_multimer_min'] = mx
                elif n == 2:
                    df.at[idx, 'ddg_multimer_min'] = round(2 * mn - mx, 6)
                else:
                    df.at[idx, 'ddg_multimer_min'] = mn
                n_fixed += 1
    if n_fixed:
        print(f"✓ Fixed {n_fixed} ddg_multimer_min NaN values")
    return df


def normalize_partners_tested(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    """Normalize FoldX partner names in partners_tested column."""
    if 'partners_tested' not in df.columns:
        return df
    n_norm = 0
    for idx in df.index:
        old = df.at[idx, 'partners_tested']
        if pd.notna(old):
            new = config.normalize_partners_tested(old, df.at[idx, 'gene'])
            if new != old:
                df.at[idx, 'partners_tested'] = new
                n_norm += 1
    if n_norm:
        print(f"✓ Normalized {n_norm} partners_tested entries")
    return df


def compute_monomer_confidence(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    """Gate monomer DDG by monomer_plddt (not best_plddt)."""
    plddt_gate = config.thresholds.get('plddt_confident', 50)

    def _conf(row):
        mp = row.get('monomer_plddt')
        dm = row.get('ddg_monomer')
        if pd.isna(mp) or pd.isna(dm):
            return False, 'low'
        if float(mp) >= plddt_gate:
            return True, 'high' if float(mp) >= 70 else 'medium'
        return False, 'low'

    results = df.apply(_conf, axis=1)
    df['ddg_monomer_confident'] = [r[0] for r in results]
    df['ddg_confidence'] = [r[1] for r in results]
    df['ddg_category'] = df['ddg_monomer'].apply(classify_ddg)
    df['ddg_category_raw'] = df['ddg_category']
    # Mark unreliable
    mask = (~df['ddg_monomer_confident']) & df['ddg_monomer'].notna()
    df.loc[mask, 'ddg_category'] = df.loc[mask, 'ddg_category'].apply(
        lambda x: f"{x}_unreliable" if x else None)
    return df


def compute_per_partner_confidence(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    """Compute confidence and classification for each per-partner DDG."""
    plddt_gate = config.thresholds.get('plddt_confident', 50)
    ddg_destab = config.thresholds.get('ddg_destabilizing', 1.0)

    for pl in config.partner_labels:
        bind_col = f"ddg_binding_{pl}"
        fold_col = f"ddg_fold_{pl}"
        plddt_col = f"multi_{pl}_plddt"
        conf_col = f"ddg_{pl}_confident"

        if bind_col not in df.columns:
            df[bind_col] = None
        if fold_col not in df.columns:
            df[fold_col] = None

        # Confidence
        df[conf_col] = df.apply(
            lambda r: pd.notna(r.get(plddt_col)) and float(sf(r.get(plddt_col))) >= plddt_gate
            if plddt_col in r.index else False, axis=1)

        # Categories
        for ddg_col, cat_col in [(bind_col, f"ddg_binding_category_{pl}"),
                                  (fold_col, f"ddg_fold_category_{pl}")]:
            if ddg_col in df.columns:
                df[cat_col] = df.apply(
                    lambda r, dc=ddg_col, cc=conf_col: (
                        classify_ddg(r[dc]) if r.get(cc, False)
                        else (f"{classify_ddg(r[dc])}_unreliable" if pd.notna(r.get(dc)) and classify_ddg(r[dc]) else None)
                    ), axis=1)

        # Interpretation
        interp_col = f"ddg_interp_{pl}"
        df[interp_col] = df.apply(
            lambda r, bc=bind_col, fc=fold_col, cc=conf_col: _interp(r, bc, fc, cc, ddg_destab),
            axis=1)

    return df


def _interp(row, bind_col, fold_col, conf_col, thresh):
    bv = row.get(bind_col)
    fv = row.get(fold_col)
    conf = row.get(conf_col, False)
    if pd.isna(bv) and pd.isna(fv):
        return None
    parts = []
    if pd.notna(fv):
        fv = float(fv)
        fd = "fold_disrupted" if fv > thresh else ("fold_stabilized" if fv < -thresh else "fold_neutral")
        parts.append(fd)
    if pd.notna(bv):
        bv = float(bv)
        bd = "binding_disrupted" if bv > thresh else ("binding_strengthened" if bv < -thresh else "binding_neutral")
        parts.append(bd)
    if not conf:
        parts.append("(low_confidence)")
    return '; '.join(parts)


def compute_ddg_summary(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    """Compute variant-level DDG summary flags and affected partners."""
    partner_labels = config.partner_labels
    ddg_destab = config.thresholds.get('ddg_destabilizing', 1.0)

    def _summary(row):
        mono = row.get('ddg_monomer')
        mono_conf = row.get('ddg_monomer_confident', False)

        if pd.notna(mono) and mono_conf:
            mv = float(mono)
            ms = "fold_disrupted" if mv > ddg_destab else ("fold_stabilized" if mv < -ddg_destab else "fold_neutral")
        elif pd.notna(mono):
            ms = "fold_unknown(low_confidence)"
        else:
            ms = "fold_unknown"

        flags = [ms]
        affected, neutral_p = [], []

        for pl in partner_labels:
            bv = row.get(f"ddg_binding_{pl}")
            conf = row.get(f"ddg_{pl}_confident", False)
            if pd.isna(bv) or not conf:
                continue
            bv = float(bv)
            if abs(bv) > ddg_destab:
                affected.append(f"{pl}({bv:+.1f})")
            else:
                neutral_p.append(pl)

        rescued, worsened = [], []
        for pl in partner_labels:
            fv = row.get(f"ddg_fold_{pl}")
            conf = row.get(f"ddg_{pl}_confident", False)
            if pd.isna(fv) or not conf:
                continue
            fv = float(fv)
            if pd.notna(mono) and mono_conf:
                mv = float(mono)
                md = "disrupted" if mv > ddg_destab else ("stabilized" if mv < -ddg_destab else "neutral")
                fd = "disrupted" if fv > ddg_destab else ("stabilized" if fv < -ddg_destab else "neutral")
                if md == "disrupted" and fd in ("neutral", "stabilized"):
                    rescued.append(pl)
                elif md in ("neutral", "stabilized") and fd == "disrupted":
                    worsened.append(pl)

        has_str = any("(-" in a for a in affected)
        has_dis = any("(+" in a for a in affected)
        if has_dis:
            flags.append("has_interaction_disrupted")
        if has_str:
            flags.append("has_interaction_strengthened")
        if has_dis and has_str:
            flags.append("mixed_interaction_effects")
        if rescued:
            flags.append(f"fold_context_rescued({';'.join(rescued)})")
        if worsened:
            flags.append(f"fold_context_worsened({';'.join(worsened)})")

        return (
            '; '.join(flags),
            ', '.join(affected) if affected else '',
            ', '.join(neutral_p) if neutral_p else '',
        )

    results = df.apply(_summary, axis=1)
    df['ddg_summary_flags'] = [r[0] for r in results]
    df['ddg_summary_partners_affected'] = [r[1] for r in results]
    df['ddg_summary_partners_neutral'] = [r[2] for r in results]
    return df


def compute_multimer_category(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    """Classify overall multimer DDG category."""
    partner_labels = config.partner_labels
    plddt_gate = config.thresholds.get('plddt_confident', 50)

    def _classify(row):
        confident_vals, all_vals = [], []
        for pl in partner_labels:
            v = row.get(f"ddg_binding_{pl}")
            if pd.notna(v):
                all_vals.append(float(v))
                if row.get(f"ddg_{pl}_confident", False):
                    confident_vals.append(float(v))
        if confident_vals:
            extreme = max(confident_vals, key=abs)
            c = classify_ddg(extreme)
            return c, c
        if all_vals:
            extreme = max(all_vals, key=abs)
            raw = classify_ddg(extreme)
            return (f"{raw}_unreliable" if raw else None), raw
        mx, mn = row.get('ddg_multimer_max'), row.get('ddg_multimer_min')
        agg = []
        if pd.notna(mx): agg.append(float(mx))
        if pd.notna(mn): agg.append(float(mn))
        if not agg:
            return None, None
        extreme = max(agg, key=abs)
        raw = classify_ddg(extreme)
        pm = row.get('multimer_plddt_max')
        if pd.notna(pm) and float(pm) < plddt_gate:
            return (f"{raw}_unreliable" if raw else None), raw
        return raw, raw

    results = df.apply(_classify, axis=1)
    df['ddg_category_multimer'] = [r[0] for r in results]
    df['ddg_category_multimer_raw'] = [r[1] for r in results]
    return df


def compute_low_confidence_flags(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    """Flag high-DDG values that lack pLDDT confidence."""
    partner_labels = config.partner_labels
    ddg_highly = config.thresholds.get('ddg_highly_destabilizing', 2.0)

    def _flag(row):
        flags = []
        for pl in partner_labels:
            bv = row.get(f"ddg_binding_{pl}")
            if pd.isna(bv):
                continue
            conf = row.get(f"ddg_{pl}_confident", False)
            if not conf and abs(float(bv)) > ddg_highly:
                p = row.get(f"multi_{pl}_plddt")
                ps = f"{int(p)}" if pd.notna(p) else "?"
                flags.append(f"{pl}:{float(bv):.2f}(pLDDT={ps})")
        return ';'.join(flags) if flags else ''

    df['ddg_low_confidence_flag'] = df.apply(_flag, axis=1)
    return df


def run_missing_monomer_ddg(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    """Run FoldX BuildModel for variants missing monomer DDG."""
    if not config.foldx_settings.get('run_missing_monomer', True):
        return df
    if config.foldx_binary is None or not config.foldx_binary.exists():
        return df

    missing = df[df['ddg_monomer'].isna()]
    if len(missing) == 0:
        return df

    print(f"Running FoldX monomer DDG for {len(missing)} missing variants...")
    foldx_dir = config.working_dir / 'foldx_expanded' / 'monomer'

    for idx, row in missing.iterrows():
        gene = str(row['gene']).lower()
        pos = int(row['position'])
        ref, alt = str(row['ref_aa']), str(row['alt_aa'])

        files = config.monomer_files(gene)
        pdb_path = files['pdb']
        if pdb_path is None:
            print(f"  ⚠ No monomer PDB for {gene} — skipping {ref}{pos}{alt}")
            continue

        work_dir = foldx_dir / f"{gene}_{ref}{pos}{alt}"
        ddg = run_buildmodel(config, pdb_path, 'A', ref, pos, alt, work_dir)
        if ddg is not None:
            df.at[idx, 'ddg_monomer'] = ddg
            df.at[idx, 'ddg_category'] = classify_ddg(ddg)
            print(f"  {gene} {ref}{pos}{alt}: DDG = {ddg:.2f}")

    print(f"  Monomer DDG coverage: {df['ddg_monomer'].notna().sum()}/{len(df)}")
    return df


def run_missing_multimer_ddg(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    """Run FoldX for interface variants missing multimer DDG."""
    if not config.foldx_settings.get('run_missing_multimer', True):
        return df
    if config.foldx_binary is None or not config.foldx_binary.exists():
        return df

    if 'is_interface_any' not in df.columns:
        return df

    missing = df[df['is_interface_any'] & df['ddg_multimer_max'].isna()]
    if len(missing) == 0:
        print("  No interface variants missing multimer DDG")
        return df

    print(f"Running FoldX multimer DDG for {len(missing)} interface variants...")
    foldx_dir = config.working_dir / 'foldx_expanded' / 'multimer'

    for idx, row in missing.iterrows():
        gene = str(row['gene']).lower()
        pos = int(row['position'])
        ref, alt = str(row['ref_aa']), str(row['alt_aa'])
        partners = str(row.get('interface_partners', ''))
        if not partners or partners == 'nan':
            continue

        for partner in partners.split(';'):
            partner = partner.strip()
            if not partner:
                continue

            # Find matching multimer definition
            matched = None
            for mdef in config.multimer_complexes:
                if mdef['gene'].lower() == gene and mdef['partner'].lower() == partner:
                    matched = (mdef['pdb'], mdef['chain_gene'], mdef['chain_partner'])
                    break
                pgm = config.partner_gene_map
                pg = pgm.get(mdef['partner'].lower(), mdef['partner'].lower())
                if mdef['gene'].lower() != gene and pg == gene and mdef['partner'].lower() == partner:
                    matched = (mdef['pdb'], mdef['chain_partner'], mdef['chain_gene'])
                    break

            if matched is None:
                continue

            pdb_file, my_chain, partner_chain = matched
            pdb_path = config.find_file(pdb_file)
            if pdb_path is None:
                continue

            query_pos = pos
            corrected, valid = config.apply_offset(gene, partner, pos)
            if not valid:
                continue
            query_pos = corrected

            work_dir = foldx_dir / f"{gene}_{ref}{pos}{alt}_{partner}"
            ddg_fold, ddg_binding = run_multimer_ddg(
                config, pdb_path, my_chain, partner_chain,
                ref, query_pos, alt, work_dir)

            if ddg_binding is not None or ddg_fold is not None:
                if ddg_binding is not None:
                    df.at[idx, f"ddg_binding_{partner}"] = ddg_binding
                if ddg_fold is not None:
                    df.at[idx, f"ddg_fold_{partner}"] = ddg_fold
                print(f"  {gene} {ref}{pos}{alt} × {partner}: "
                      f"fold={ddg_fold}, binding={ddg_binding}")

    return df
