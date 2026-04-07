"""Disruption scoring, tier assignment, mechanism classification, and concordance."""

import pandas as pd
import numpy as np

from .config import PipelineConfig
from .constants import (
    classify_ddg, sf, ss, sb,
    BURIAL_RANK, RANK_TO_BURIAL,
)


# ─── Franklin standardization ────────────────────────────────────────────────

def std_franklin(v):
    if pd.isna(v):
        return 'No data'
    v = str(v).strip()
    vl = v.lower()
    if 'pathogenic' in vl and 'likely' in vl:
        return 'Likely Pathogenic'
    elif 'pathogenic' in vl:
        return 'Pathogenic'
    elif 'benign' in vl and 'likely' in vl:
        return 'Likely Benign'
    elif 'benign' in vl:
        return 'Benign'
    elif 'vus' in vl and 'high' in vl:
        return 'VUS (high)'
    elif 'vus' in vl and 'mid' in vl:
        return 'VUS (mid)'
    elif 'vus' in vl and 'low' in vl:
        return 'VUS (low)'
    elif 'vus' in vl:
        return 'VUS'
    return v


def classify_am(val, thresholds=None):
    if thresholds is None:
        thresholds = {'pathogenic': 0.564, 'benign': 0.340}
    if pd.isna(val):
        return val
    s = str(val).strip()
    if s in ('likely_pathogenic', 'likely_benign', 'ambiguous'):
        return s
    try:
        score = float(s)
        if score >= thresholds['pathogenic']:
            return 'likely_pathogenic'
        elif score < thresholds['benign']:
            return 'likely_benign'
        return 'ambiguous'
    except ValueError:
        return s


# ─── Pipeline 1 Scoring ──────────────────────────────────────────────────────

def compute_score(row, config: PipelineConfig):
    """Pipeline 1 single-residue disruption score."""
    t = config.thresholds
    partner_labels = config.partner_labels
    disruption_pts = [(p[0], p[1]) for p in t['disruption_points']]
    plddt_gate = t.get('plddt_confident', 50)
    contact_thresh = t.get('contact_driven_threshold', 6)

    score, ev = 0.0, []
    sev = sf(row.get('substitution_severity'), 0)
    mono_c = sf(row.get('monomer_n_contacts'), 0)

    max_inter = 0.0
    best_inter_partner = None
    gated_iface_partners = []

    for pl in partner_labels:
        plddt_col = f"multi_{pl}_plddt"
        ic_col = f"multi_{pl}_inter_contacts"
        iface_col = f"multi_{pl}_is_interface"

        pl_plddt = sf(row.get(plddt_col), 0) if plddt_col in row.index else 0
        if pl_plddt < plddt_gate:
            continue
        if iface_col in row.index and sb(row.get(iface_col)):
            gated_iface_partners.append(pl)
        if ic_col in row.index and pd.notna(row[ic_col]):
            ic_val = float(row[ic_col])
            if ic_val > max_inter:
                max_inter = ic_val
                best_inter_partner = pl

    total_contacts = mono_c + max_inter
    disruption = round(sev * total_contacts, 2)

    for thresh, pts in disruption_pts:
        if disruption >= thresh:
            score += pts
            ev.append(f'disruption({disruption:.1f})')
            break
    else:
        ev.append(f'no_disruption({disruption:.1f})')

    # Interface bonus
    n_gated = len(gated_iface_partners)
    if n_gated >= 2:
        score += 2.0
        ev.append(f'multi_interface({n_gated})')
    elif n_gated == 1:
        score += 1.5
        ev.append(f'interface({gated_iface_partners[0]})')

    # Burial
    mono_burial = ss(row.get('monomer_burial'))
    best_rank = BURIAL_RANK.get(mono_burial, 0)
    best_source = 'monomer'
    for pl in partner_labels:
        burial_col = f"multi_{pl}_burial"
        plddt_col = f"multi_{pl}_plddt"
        if burial_col in row.index and pd.notna(row.get(burial_col)):
            if sf(row.get(plddt_col), 0) >= plddt_gate:
                pl_rank = BURIAL_RANK.get(ss(row[burial_col]), 0)
                if pl_rank > best_rank:
                    best_rank = pl_rank
                    best_source = pl
    best_burial = RANK_TO_BURIAL.get(best_rank, 'unknown')
    if best_burial == 'buried_core':
        score += 2.0
        ev.append(f'buried_core({best_source})')
    elif best_burial == 'partially_buried':
        score += 1.0
        ev.append(f'partially_buried({best_source})')

    # pLDDT multiplier
    bp = row.get('best_plddt')
    if bp is not None and not pd.isna(bp):
        bp_val = float(bp)
        if bp_val < 50:
            score *= 0.4
            ev.append(f'very_low_plddt({int(bp_val)})')
        elif bp_val < 70:
            score *= 0.7
            ev.append(f'low_plddt({int(bp_val)})')

    return pd.Series({
        'v6_final_score': round(score, 2),
        'v6_score_evidence': ';'.join(ev),
        'v6_contact_disruption': disruption,
        'v6_total_contacts': total_contacts,
        'v6_max_inter_contacts': max_inter,
        'v6_best_inter_partner': best_inter_partner or '',
        'v6_best_burial': best_burial,
        'v6_best_burial_source': best_source,
        'v6_interface_partners_gated': ';'.join(gated_iface_partners) if gated_iface_partners else '',
        'v6_n_interface_partners_gated': n_gated,
    })


def assign_tier(score, config: PipelineConfig):
    t = config.thresholds
    if pd.isna(score):
        return t.get('default_tier', 'Tier 4')
    s = float(score)
    for thresh, label in t['tiers']:
        if s >= thresh:
            return label
    return t.get('default_tier', 'Tier 4')


# ─── Evaluability ────────────────────────────────────────────────────────────

def compute_evaluability(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    partner_labels = config.partner_labels
    plddt_gate = config.thresholds.get('plddt_confident', 50)

    df['structure_evaluable'] = df['best_plddt'].apply(
        lambda x: pd.notna(x) and float(x) >= plddt_gate)
    df['ddg_evaluable'] = df['ddg_confidence'].apply(
        lambda x: ss(x).lower() not in ('low', ''))
    df['am_evaluable'] = df['AlphaMissense'].notna() if 'AlphaMissense' in df.columns else False
    df['franklin_evaluable'] = df['franklin'].notna() if 'franklin' in df.columns else False

    def _check_multi(row):
        for pl in partner_labels:
            col = f"multi_{pl}_plddt"
            if col in row.index and pd.notna(row[col]) and float(row[col]) >= plddt_gate:
                return True
        return False

    df['multimer_evaluable'] = df.apply(_check_multi, axis=1)
    return df


# ─── 16-Category Mechanism Classification ────────────────────────────────────

def classify_mechanism(row, config: PipelineConfig):
    """16-category mechanism classification."""
    partner_labels = config.partner_labels
    ddg_destab = config.thresholds.get('ddg_destabilizing', 1.0)
    contact_thresh = config.thresholds.get('contact_driven_threshold', 6)

    tier = ss(row.get('v6_tier'))
    is_high_tier = tier in ('Tier 1', 'Tier 2')
    structure_eval = row.get('structure_evaluable', False)

    ddg_m = row.get('ddg_monomer')
    mono_conf = row.get('ddg_monomer_confident', False)
    mono_c = sf(row.get('monomer_n_contacts'), 0)
    gated_iface = ss(row.get('v6_interface_partners_gated'))
    is_iface = len(gated_iface) > 0

    # Fold axis
    has_dm = pd.notna(ddg_m)
    dm_val = float(ddg_m) if has_dm else 0.0

    mono_fold_destab = has_dm and mono_conf and dm_val > ddg_destab
    mono_fold_stab = has_dm and mono_conf and dm_val < -ddg_destab

    partner_fold_destab = False
    partner_fold_stab = False
    for pl in partner_labels:
        fv = row.get(f"ddg_fold_{pl}")
        conf = row.get(f"ddg_{pl}_confident", False)
        if pd.notna(fv) and conf:
            if float(fv) > ddg_destab:
                partner_fold_destab = True
            if float(fv) < -ddg_destab:
                partner_fold_stab = True

    fold_destab = mono_fold_destab or partner_fold_destab
    fold_stab = mono_fold_stab or partner_fold_stab

    # PPI axis
    ppi_destab, ppi_stab = False, False
    ppi_partner_destab, ppi_partner_stab = None, None

    has_per_partner = False
    for pl in partner_labels:
        bv = row.get(f"ddg_binding_{pl}")
        if pd.notna(bv):
            has_per_partner = True
            conf = row.get(f"ddg_{pl}_confident", False)
            if not conf:
                continue
            bv = float(bv)
            if bv > ddg_destab:
                ppi_destab = True
                if ppi_partner_destab is None:
                    ppi_partner_destab = pl
            if bv < -ddg_destab:
                ppi_stab = True
                if ppi_partner_stab is None:
                    ppi_partner_stab = pl

    # Fallback to aggregates
    if not has_per_partner:
        ddg_x = row.get('ddg_multimer_max')
        ddg_n = row.get('ddg_multimer_min')
        pt = ss(row.get('partners_tested'))
        tested = set(p.strip() for p in pt.split(';') if p.strip()) if pt else set()
        iface_set = set(gated_iface.split(';')) if gated_iface else set()

        if tested:
            dx = pd.notna(ddg_x) and float(ddg_x) > ddg_destab
            dn = pd.notna(ddg_n) and float(ddg_n) > ddg_destab
            sx = pd.notna(ddg_x) and float(ddg_x) < -ddg_destab
            sn = pd.notna(ddg_n) and float(ddg_n) < -ddg_destab
            overlap = iface_set & tested if is_iface else set()

            if (dx or dn) and (overlap or is_iface):
                ppi_destab = True
                ppi_partner_destab = sorted(overlap)[0] if overlap else (sorted(iface_set)[0] if iface_set else '')
            if (sx or sn) and (overlap or is_iface):
                ppi_stab = True
                ppi_partner_stab = sorted(overlap)[0] if overlap else (sorted(iface_set)[0] if iface_set else '')

    ppi_partner = ppi_partner_destab or ppi_partner_stab or ''

    # Structure unevaluable
    if not structure_eval:
        am = ss(row.get('AlphaMissense')).lower()
        fr = std_franklin(row.get('franklin')).lower()
        ext = (am == 'likely_pathogenic') or (fr in ('pathogenic', 'likely pathogenic', 'vus (high)'))
        return 'Structure unevaluable', ppi_partner, ext

    # Fold disrupted
    if fold_destab:
        if ppi_destab:
            return 'Fold + PPI destabilization', ppi_partner_destab or '', False
        elif ppi_stab:
            return 'Fold destab. + PPI stabilization (conflicting)', ppi_partner_stab or '', False
        elif is_iface:
            return 'Fold destab. at interface', '', False
        return 'Fold destabilization', '', False

    # Fold stabilized
    if fold_stab:
        if ppi_stab:
            return 'Stabilizing + PPI (potential GoF)', ppi_partner_stab or '', False
        elif ppi_destab:
            return 'Stabilizing + PPI destabilization (conflicting)', ppi_partner_destab or '', False
        elif is_iface:
            return 'Stabilizing at interface (potential GoF)', '', False
        return 'Stabilizing (potential GoF)', '', False

    # Fold neutral
    if ppi_destab and not ppi_stab:
        return 'PPI destabilization', ppi_partner_destab or '', False
    if ppi_stab and not ppi_destab:
        return 'PPI stabilization (potential GoF)', ppi_partner_stab or '', False
    if ppi_destab and ppi_stab:
        return 'PPI conflicting (mixed partner signals)', ppi_partner, False
    if is_iface:
        return 'Interface variant (DDG neutral)', '', False
    if is_high_tier:
        if mono_c >= contact_thresh:
            return 'Structural variant - contact-driven (DDG neutral)', '', False
        return 'Structural variant - burial-driven (DDG neutral)', '', False

    return 'Benign (structurally evaluated)', '', False


# ─── Four-Way Concordance ────────────────────────────────────────────────────

def compute_concordance(row, config: PipelineConfig):
    """Four-way concordance: strict, relaxed, T3-inclusive."""
    t = config.thresholds
    partner_labels = config.partner_labels
    ddg_destab = t.get('ddg_destabilizing', 1.0)
    ddg_highly = t.get('ddg_highly_destabilizing', 2.0)
    am_path_thresh = t.get('am_pathogenic', 0.564)

    tier = ss(row.get('v6_tier'))
    am = ss(row.get('AlphaMissense')).lower()
    fr = std_franklin(row.get('franklin'))
    fr_lower = fr.lower() if isinstance(fr, str) else ''
    ddg_conf = ss(row.get('ddg_confidence')).lower()

    # Three-axis max |DDG|
    ddg_vals = []
    if pd.notna(row.get('ddg_monomer')) and row.get('ddg_monomer_confident', False):
        ddg_vals.append(abs(float(row['ddg_monomer'])))

    has_pp = False
    for pl in partner_labels:
        for col_pfx in ['ddg_binding', 'ddg_fold']:
            col = f"{col_pfx}_{pl}"
            if col in row.index and pd.notna(row.get(col)):
                has_pp = True
                if row.get(f"ddg_{pl}_confident", False):
                    ddg_vals.append(abs(float(row[col])))

    if not has_pp:
        for col in ['ddg_multimer_max', 'ddg_multimer_min']:
            v = row.get(col)
            if pd.notna(v):
                ddg_vals.append(abs(float(v)))

    max_abs_ddg = max(ddg_vals) if ddg_vals else 0.0

    struct_eval = row.get('structure_evaluable', False)
    ddg_eval = row.get('ddg_evaluable', False)
    am_eval = row.get('am_evaluable', False)
    fr_eval = row.get('franklin_evaluable', False)

    tier_t12 = tier in ('Tier 1', 'Tier 2')
    tier_t123 = tier_t12 or tier == 'Tier 3'

    ddg_strict = 1 if (ddg_conf == 'high' and max_abs_ddg >= ddg_highly) else 0
    ddg_relaxed = 1 if (ddg_conf not in ('low', '') and max_abs_ddg >= ddg_destab) else 0
    am_strict = 1 if am == 'likely_pathogenic' else 0
    am_relaxed = 1 if am in ('likely_pathogenic', 'ambiguous') else 0
    fr_strict = 1 if fr_lower in ('pathogenic', 'likely pathogenic', 'vus (high)') else 0
    fr_relaxed = 1 if fr_lower in ('pathogenic', 'likely pathogenic', 'vus (high)', 'vus (mid)') else 0

    def build(tier_v, ddg_v, am_v, fr_v):
        s, d = 0, 0
        if struct_eval:
            s += tier_v; d += 1
        if ddg_eval:
            s += ddg_v; d += 1
        if am_eval:
            s += am_v; d += 1
        if fr_eval:
            s += fr_v; d += 1
        return s, max(d, 1)

    s_s, s_d = build(1 if tier_t12 else 0, ddg_strict, am_strict, fr_strict)
    r_s, r_d = build(1 if tier_t12 else 0, ddg_relaxed, am_relaxed, fr_relaxed)
    t3_s, t3_d = build(1 if tier_t123 else 0, ddg_strict, am_strict, fr_strict)

    return pd.Series({
        'four_way_strict': s_s, 'four_way_strict_denom': s_d,
        'concordance_strict': f"{s_s}/{s_d}",
        'four_way_relaxed': r_s, 'four_way_relaxed_denom': r_d,
        'concordance_relaxed': f"{r_s}/{r_d}",
        'four_way_t3': t3_s, 'four_way_t3_denom': t3_d,
        'concordance_t3': f"{t3_s}/{t3_d}",
        'structure_vote_strict': 1 if tier_t12 else 0,
        'structure_vote_t3': 1 if tier_t123 else 0,
        'ddg_vote_strict': ddg_strict, 'ddg_vote_relaxed': ddg_relaxed,
        'am_vote_strict': am_strict, 'am_vote_relaxed': am_relaxed,
        'franklin_vote_strict': fr_strict, 'franklin_vote_relaxed': fr_relaxed,
        'max_abs_ddg': max_abs_ddg,
    })


# ─── Annotation normalization ────────────────────────────────────────────────

def normalize_annotations(df: pd.DataFrame, config: PipelineConfig) -> pd.DataFrame:
    """Normalize AlphaMissense and Franklin annotations."""
    t = config.thresholds
    am_thresh = {'pathogenic': t.get('am_pathogenic', 0.564),
                 'benign': t.get('am_benign', 0.340)}

    if 'AlphaMissense' in df.columns:
        if 'AlphaMissense_raw' not in df.columns:
            df['AlphaMissense_raw'] = df['AlphaMissense']
        df['AlphaMissense'] = df['AlphaMissense'].apply(lambda v: classify_am(v, am_thresh))

    if 'franklin' in df.columns:
        if 'franklin_raw' not in df.columns:
            df['franklin_raw'] = df['franklin']
        df['franklin'] = df['franklin'].apply(
            lambda v: v if pd.isna(v) else std_franklin(v))

    return df
