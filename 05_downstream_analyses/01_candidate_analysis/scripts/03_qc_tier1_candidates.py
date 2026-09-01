#!/usr/bin/env python3

# ============================================================
# QC Tier-1 CRE candidates
#
# Purpose:
#   Validate focal and secondary Tier-1 candidates against
#   detailed per-species CRE-classification results.
#
# Checks:
#   - expected present/turnover_candidate state pattern
#   - successful homologous-locus mapping
#   - reciprocal positional overlap for present CREs
#   - local same-FBgn support for turnover candidates
#   - consistency between focal and singleton analyses
#
# Input/output paths and QC parameters:
#   Supplied by the pipeline wrapper using
#   config/candidate_config.sh.
# ============================================================

import argparse
from datetime import datetime
from pathlib import Path
import platform
import sys

import pandas as pd


# ============================================================
# Arguments
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description='QC focal and secondary Tier-1 CRE candidates against detailed per-species CRE classifications.')
    parser.add_argument('--focal-dir', type=Path, required=True)
    parser.add_argument('--secondary', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--turnover-dir', type=Path, required=True)
    parser.add_argument('--out-summary', type=Path, required=True)
    parser.add_argument('--out-detail', type=Path, required=True)
    parser.add_argument('--metadata-out', type=Path, required=True)
    parser.add_argument('--min-reciprocal-overlap', type=float, default=0.5)
    parser.add_argument('--expected-reference-cres', type=int, default=None)
    parser.add_argument('--recurrence-min-clades', type=int, default=2)
    return parser.parse_args()

# ============================================================
# CRE-state definitions
# ============================================================

VALID_TIER1_STATES = {"present", "turnover_candidate"}

# ============================================================
# Required columns
# ============================================================

FOCAL_REQUIRED_COLUMNS = {
    "group_name",
    "dmel_cre_id",
    "category",
    "focal_species",
    "focal_state",
    "comparison_species",
    "comparison_consensus_state",
}

SECONDARY_REQUIRED_COLUMNS = {
    "dmel_cre_id",
    "group_name",
    "discordant_species",
    "discordant_state",
    "consensus_state",
    "group_species",
}

TURNOVER_REQUIRED_COLUMNS = {
    "dmel_cre_id",
    "class",
    "alignment_status",
    "target_chrom",
    "target_start0",
    "target_end0",
    "best_peak_id",
    "best_overlap_bp",
    "best_overlap_fraction_lifted",
    "best_overlap_fraction_peak",
    "positional_peak",
    "same_fbgn_peak_local",
    "local_same_fbgn_peak_id",
    "local_same_fbgn",
    "best_local_same_fbgn_distance_bp",
    "shared_fbgn_with_best_peak",
    "gene_support_at_best_peak",
}

# ============================================================
# Helpers
# ============================================================

def require_file(path):
    if not path.is_file():
        raise SystemExit(f'ERROR: required file not found:\n{path}')

def require_directory(path):
    if not path.is_dir():
        raise SystemExit(f'ERROR: required directory not found:\n{path}')

def require_columns(df, required, label):
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise SystemExit(f'ERROR: missing columns in {label}:\n' + '\n'.join(missing))

def split_pipe(value):
    return [x for x in str(value).split('|') if x]

def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0

def is_missing(value):
    if pd.isna(value):
        return True
    return str(value).strip().lower() in {'', 'na', 'nan', 'none'}

def validate_species_row(row, expected_state, min_reciprocal_overlap):
    """
    Validate one CRE x species classification record.

    Tier-1 allows only:
        present
        turnover_candidate
    """
    flags = []
    observed_class = str(row.get('class', '')).strip()
    alignment_status = str(row.get('alignment_status', '')).strip()
    mapped = alignment_status == 'mapped'
    frac_lifted = as_float(row.get('best_overlap_fraction_lifted', 0))
    frac_peak = as_float(row.get('best_overlap_fraction_peak', 0))
    positional = str(row.get('positional_peak', '')).strip() == 'yes'
    same_fbgn_local = str(row.get('same_fbgn_peak_local', '')).strip() == 'yes'
    local_peak_id = row.get('local_same_fbgn_peak_id', '')

    if observed_class != expected_state:
        flags.append('STATE_MISMATCH')

    if not mapped:
        flags.append('NOT_MAPPED')
        
    if observed_class == 'present':
        if not positional:
            flags.append('PRESENT_WITHOUT_POSITIONAL_PEAK')
        if frac_lifted < min_reciprocal_overlap:
            flags.append('PRESENT_LIFTED_OVERLAP_BELOW_THRESHOLD')
        if frac_peak < min_reciprocal_overlap:
            flags.append('PRESENT_PEAK_OVERLAP_BELOW_THRESHOLD')

    elif observed_class == 'turnover_candidate':
        if positional:
            flags.append('TURNOVER_HAS_POSITIONAL_PEAK')
        if not same_fbgn_local:
            flags.append('TURNOVER_WITHOUT_LOCAL_SAME_FBGN')
        if is_missing(local_peak_id):
            flags.append('TURNOVER_WITHOUT_LOCAL_PEAK_ID')

    else:
        flags.append(f'UNEXPECTED_CLASS_{observed_class}')
    return sorted(set(flags))


# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()
    if not 0.0 <= args.min_reciprocal_overlap <= 1.0:
        raise SystemExit('ERROR: min-reciprocal-overlap must be between 0 and 1.')
    if args.recurrence_min_clades < 1:
        raise SystemExit('ERROR: recurrence-min-clades must be >= 1.')
    require_directory(args.focal_dir)
    require_file(args.secondary)
    require_file(args.reference)
    require_directory(args.turnover_dir)
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    args.out_detail.parent.mkdir(parents=True, exist_ok=True)
    args.metadata_out.parent.mkdir(parents=True, exist_ok=True)

    # ========================================================
    # Load focal Tier-1 candidates
    # ========================================================

    focal_files = sorted(args.focal_dir.glob('*_tier1.tsv'))
    if not focal_files:
        raise SystemExit(f'ERROR: no focal *_tier1.tsv files found in:\n{args.focal_dir}')
    focal_tables = []
    for path in focal_files:
        df = pd.read_csv(path, sep='\t', dtype=str).fillna('')
        require_columns(df, FOCAL_REQUIRED_COLUMNS, str(path))
        if df.empty:
            continue
        if not (df['category'] == 'tier1').all():
            raise SystemExit(f'ERROR: non-tier1 rows found in:\n{path}')
        focal_tables.append(df)
    if focal_tables:
        focal = pd.concat(focal_tables, ignore_index=True)
    else:
        focal = pd.DataFrame(columns=sorted(FOCAL_REQUIRED_COLUMNS))
        
    # --------------------------------------------------------
    # Validate focal uniqueness
    # --------------------------------------------------------

    if not focal.empty:
        if focal.duplicated(subset=['group_name', 'dmel_cre_id']).any():
            raise SystemExit('ERROR: duplicate focal group x CRE records.')

    # ========================================================
    # Load secondary/singleton Tier-1 candidates
    # ========================================================

    secondary = pd.read_csv(args.secondary, sep='\t', dtype=str).fillna('')
    require_columns(secondary, SECONDARY_REQUIRED_COLUMNS, 'secondary Tier-1 table')
    if not secondary.empty:
        if secondary.duplicated(subset=['group_name', 'dmel_cre_id']).any():
            raise SystemExit('ERROR: duplicate secondary group x CRE records.')

    # ========================================================
    # Convert focal records into common evidence representation
    # ========================================================

    evidence = {}
    for _, row in focal.iterrows():
        group = row['group_name']
        cre = row['dmel_cre_id']
        focal_species = row['focal_species']
        comparisons = split_pipe(row['comparison_species'])
        focal_state = row['focal_state']
        consensus_state = row['comparison_consensus_state']
        expected_states = {focal_species: focal_state}
        expected_states.update({sp: consensus_state for sp in comparisons})
        evidence[group, cre] = {'group_name': group, 
                                'dmel_cre_id': cre, 
                                'species': [focal_species, *comparisons], 
                                'expected_states': expected_states, 
                                'is_focal_tier1': True, 
                                'is_singleton_tier1': False, 
                                'focal_species': focal_species, 
                                'focal_state': focal_state,
                                'comparison_species': '|'.join(comparisons), 
                                'comparison_consensus_state': consensus_state, 
                                'discordant_species': focal_species, 
                                'discordant_state': focal_state, 
                                'consensus_state': consensus_state}

    # ========================================================
    # Add singleton Tier-1 evidence
    # ========================================================

    for _, row in secondary.iterrows():
        group = row['group_name']
        cre = row['dmel_cre_id']
        group_species = split_pipe(row['group_species'])
        discordant_species = row['discordant_species']
        discordant_state = row['discordant_state']
        consensus_state = row['consensus_state']
        expected_states = {sp: discordant_state if sp == discordant_species else consensus_state for sp in group_species}
        key = (group, cre)

        # ----------------------------------------------------
        # Same group x CRE may already exist as focal Tier-1.
        # In that case, both analyses must imply exactly
        # the same state pattern.
        # ----------------------------------------------------

        if key in evidence:
            existing = evidence[key]
            if existing['expected_states'] != expected_states:
                raise SystemExit(f"ERROR: focal and singleton Tier-1 state patterns disagree for:\n{group} / {cre}\n\nFocal: {existing['expected_states']}\nSingleton: {expected_states}")
            existing['is_singleton_tier1'] = True
            existing['discordant_species'] = discordant_species
            existing['discordant_state'] = discordant_state
            existing['consensus_state'] = consensus_state
        else:
            evidence[key] = {'group_name': group, 
                             'dmel_cre_id': cre, 
                             'species': group_species, 
                             'expected_states': expected_states, 
                             'is_focal_tier1': False, 
                             'is_singleton_tier1': True, 
                             'focal_species': 'NA', 
                             'focal_state': 'NA', 
                             'comparison_species': 'NA', 
                             'comparison_consensus_state': 'NA', 
                             'discordant_species': discordant_species, 
                             'discordant_state': discordant_state, 
                             'consensus_state': consensus_state}
    # ========================================================
    # Cross-step consistency:
    # every focal Tier-1 should also be a singleton Tier-1
    # ========================================================

    focal_not_singleton = [(record['group_name'], record['dmel_cre_id']) for record in evidence.values() if record['is_focal_tier1'] and (not record['is_singleton_tier1'])]
    if focal_not_singleton:
        lines = [f'{group}\t{cre}' for group, cre in focal_not_singleton]
        raise SystemExit('ERROR: focal Tier-1 candidates were not found by the singleton Tier-1 analysis.\nStep 01 and Step 02 are inconsistent:\n' + '\n'.join(lines))

    # ========================================================
    # Assign evidence type
    # ========================================================

    for record in evidence.values():
        if record['is_focal_tier1']:
            record['tier1_scope'] = 'focal'
        else:
            record['tier1_scope'] = 'secondary_only'

    # ========================================================
    # Determine all species required
    # ========================================================

    species_needed = sorted({sp for record in evidence.values() for sp in record['species']})

    # ========================================================
    # Load per-species turnover tables once
    # ========================================================

    for species in species_needed:
        path = args.turnover_dir / f'dmel_to_{species}_cre_turnover.tsv'
        require_file(path)
        df = pd.read_csv(path, sep='\t', dtype=str).fillna('')
        require_columns(df, TURNOVER_REQUIRED_COLUMNS, str(path))
        if df['dmel_cre_id'].duplicated().any():
            raise SystemExit(f'ERROR: duplicate dmel_cre_id values in:\n{path}')
        turnover[species] = df.set_index('dmel_cre_id')

    # ========================================================
    # Load Dmel reference annotation
    # ========================================================

    ref = pd.read_csv(args.reference, sep='\t', dtype=str).fillna('')
    if args.expected_reference_cres is not None and len(ref) != args.expected_reference_cres:
        raise SystemExit(f'ERROR: unexpected number of D. melanogaster reference CREs.\nExpected: {args.expected_reference_cres}\nObserved: {len(ref)}')
    if 'dmel_cre_id' not in ref.columns:
        raise SystemExit('ERROR: Dmel reference table lacks dmel_cre_id.')
    if ref['dmel_cre_id'].duplicated().any():
        raise SystemExit('ERROR: duplicate CRE IDs in Dmel reference table.')
    ref = ref.set_index('dmel_cre_id')

    # ========================================================
    # QC all unique group x CRE evidence records
    # ========================================================

    detail_rows = []
    summary_rows = []

    records = sorted(evidence.values(), key=lambda x: (x['group_name'], x['dmel_cre_id']))
    for record in records:
        group = record['group_name']
        cre = record['dmel_cre_id']
        species = record['species']
        expected_states = record['expected_states']
        candidate_flags = []
        
        # ----------------------------------------------------
        # Tier-1 pattern itself
        # ----------------------------------------------------

        unique_expected_states = set(expected_states.values())
        if not unique_expected_states.issubset(VALID_TIER1_STATES):
            candidate_flags.append('INVALID_TIER1_STATE')
        if unique_expected_states != {'present', 'turnover_candidate'}:
            candidate_flags.append('NOT_PRESENT_VS_TURNOVER')

        # ----------------------------------------------------
        # Dmel annotation
        # ----------------------------------------------------

        if cre in ref.index:
            ref_row = ref.loc[cre]
            dmel_fbgn = ref_row.get('fbgn_target_genes', 'NA')
            chrom = ref_row.get('chrom', 'NA')
            start0 = ref_row.get('start0', 'NA')
            end0 = ref_row.get('end0', 'NA')
            training_set = ref_row.get('training_set', 'NA')
            method = ref_row.get('method', 'NA')
            scrmshaw_score = ref_row.get('scrmshaw_score', 'NA')
            rank = ref_row.get('rank', 'NA')
        else:
            dmel_fbgn = 'NA'
            chrom = 'NA'
            start0 = 'NA'
            end0 = 'NA'
            training_set = 'NA'
            method = 'NA'
            scrmshaw_score = 'NA'
            rank = 'NA'
            candidate_flags.append('MISSING_DMel_REFERENCE_ANNOTATION')

        # ----------------------------------------------------
        # Species-level QC
        # ----------------------------------------------------

        n_species_pass = 0
        for sp in species:
            expected_state = expected_states[sp]
            if cre not in turnover[sp].index:
                species_flags = ['CRE_MISSING_FROM_SPECIES_TABLE']
                row = None
            else:
                row = turnover[sp].loc[cre]
                species_flags = validate_species_row(row, expected_state, args.min_reciprocal_overlap)
            if not species_flags:
                n_species_pass += 1
            else:
                candidate_flags.append(f'SPECIES_QC_FAIL_{sp}')

            # ------------------------------------------------
            # Species role
            # ------------------------------------------------

            if sp == record['discordant_species']:
                role = 'discordant'
            else:
                role = 'consensus'

            # ------------------------------------------------
            # Detail row
            # ------------------------------------------------

            detail_rows.append({'group_name': group, 
                                'dmel_cre_id': cre, 
                                'tier1_scope': record['tier1_scope'], 
                                'is_focal_tier1': 'yes' if record['is_focal_tier1'] else 'no', 
                                'is_singleton_tier1': 'yes' if record['is_singleton_tier1'] else 'no', 
                                'species': sp, 
                                'role': role, 
                                'expected_state': expected_state, 
                                'observed_class': row.get('class', 'NA') if row is not None else 'NA', 
                                'alignment_status': row.get('alignment_status', 'NA') if row is not None else 'NA', 
                                'target_chrom': row.get('target_chrom', 'NA') if row is not None else 'NA', 
                                'target_start0': row.get('target_start0', 'NA') if row is not None else 'NA', 
                                'target_end0': row.get('target_end0', 'NA') if row is not None else 'NA', 
                                'best_peak_id': row.get('best_peak_id', 'NA') if row is not None else 'NA', 
                                'best_overlap_bp': row.get('best_overlap_bp', 'NA') if row is not None else 'NA', 
                                'best_overlap_fraction_lifted': row.get('best_overlap_fraction_lifted', 'NA') if row is not None else 'NA', 
                                'best_overlap_fraction_peak': row.get('best_overlap_fraction_peak', 'NA') if row is not None else 'NA', 
                                'positional_peak': row.get('positional_peak', 'NA') if row is not None else 'NA', 
                                'same_fbgn_peak_local': row.get('same_fbgn_peak_local', 'NA') if row is not None else 'NA', 
                                'local_same_fbgn_peak_id': row.get('local_same_fbgn_peak_id', 'NA') if row is not None else 'NA', 
                                'local_same_fbgn': row.get('local_same_fbgn', 'NA') if row is not None else 'NA', 
                                'best_local_same_fbgn_distance_bp': row.get('best_local_same_fbgn_distance_bp', 'NA') if row is not None else 'NA', 
                                'shared_fbgn_with_best_peak': row.get('shared_fbgn_with_best_peak', 'NA') if row is not None else 'NA', 
                                'gene_support_at_best_peak': row.get('gene_support_at_best_peak', 'NA') if row is not None else 'NA', 
                                'species_qc': 'PASS' if not species_flags else ';'.join(species_flags)})

        # ----------------------------------------------------
        # Candidate-level QC
        # ----------------------------------------------------

        n_expected_species = len(species)
        if n_species_pass != n_expected_species:
            candidate_flags.append('SPECIES_LEVEL_QC_FAIL')
        candidate_flags = sorted(set(candidate_flags))
        final_qc = 'PASS' if not candidate_flags else ';'.join(candidate_flags)

        # ----------------------------------------------------
        # Summary row
        # ----------------------------------------------------

        summary_rows.append({
            "group_name": group,
            "dmel_cre_id": cre,
            "tier1_scope": record["tier1_scope"],
            "is_focal_tier1": "yes" if record["is_focal_tier1"] else "no",
            "is_singleton_tier1": "yes" if record["is_singleton_tier1"] else "no",
            "discordant_species": record["discordant_species"],
            "discordant_state": record["discordant_state"],
            'consensus_state': record['consensus_state'], 
            'group_species': '|'.join(species), 
            'n_group_species': len(species), 
            'dmel_fbgn_target_genes': dmel_fbgn, 
            'dmel_chrom': chrom, 
            'dmel_start0': start0, 
            'dmel_end0': end0, 
            'training_set': training_set, 
            'method': method, 
            'dmel_scrmshaw_score': scrmshaw_score, 
            'dmel_rank': rank, 
            'n_species_expected': n_expected_species, 
            'n_species_passing_qc': n_species_pass, 
            'candidate_qc': final_qc})

    # ========================================================
    # Build output tables
    # ========================================================

    summary = pd.DataFrame(summary_rows)
    detail = pd.DataFrame(detail_rows)

    # --------------------------------------------------------
    # Stable sorting
    # --------------------------------------------------------

    if not summary.empty:
        summary = summary.sort_values(['tier1_scope', 'group_name', 'dmel_cre_id'], kind='mergesort').reset_index(drop=True)
    if not detail.empty:
        detail = detail.sort_values(['tier1_scope', 'group_name', 'dmel_cre_id', 'role', 'species'], kind='mergesort').reset_index(drop=True)

    # --------------------------------------------------------
    # Final uniqueness checks
    # --------------------------------------------------------

    if not summary.empty:
        if summary.duplicated(subset=['group_name', 'dmel_cre_id']).any():
            raise SystemExit('ERROR: duplicate group x CRE rows in QC summary.')
    if not detail.empty:
        if detail.duplicated(subset=['group_name', 'dmel_cre_id', 'species']).any():
            raise SystemExit('ERROR: duplicate group x CRE x species rows in QC detail table.')

    # ========================================================
    # Write outputs
    # ========================================================

    summary.to_csv(args.out_summary, sep='\t', index=False)
    detail.to_csv(args.out_detail, sep='\t', index=False)
    
    # ========================================================
    # Metadata
    # ========================================================

    metadata = pd.DataFrame([{'script': Path(__file__).name, 
                              'run_timestamp': datetime.now().astimezone().isoformat(), 
                              'python_version': sys.version.split()[0], 
                              'pandas_version': pd.__version__, 
                              'platform': platform.platform(), 
                              'min_reciprocal_overlap': args.min_reciprocal_overlap, 
                              'recurrence_min_clades': args.recurrence_min_clades, 
                              'expected_reference_cres': args.expected_reference_cres if args.expected_reference_cres is not None else 'NA', 
                              'n_focal_tier1_rows': len(focal), 
                              'n_singleton_tier1_rows': len(secondary), 
                              'n_unique_group_cre_candidates': len(summary), 
                              'n_unique_candidate_cres': summary['dmel_cre_id'].nunique() if not summary.empty else 0, 
                              'n_focal_candidates': (summary['tier1_scope'] == 'focal').sum() if not summary.empty else 0, 
                              'n_secondary_only_candidates': (summary['tier1_scope'] == 'secondary_only').sum() if not summary.empty else 0, 
                              'n_qc_pass': (summary['candidate_qc'] == 'PASS').sum() if not summary.empty else 0, 
                              'focal_input_dir': str(args.focal_dir.resolve()), 
                              'secondary_input': str(args.secondary.resolve()), 
                              'reference_input': str(args.reference.resolve()), 
                              'turnover_dir': str(args.turnover_dir.resolve())}])
    metadata.to_csv(args.metadata_out, sep='\t', index=False)

    # ========================================================
    # Console summary
    # ========================================================

     print()
    print('=' * 72)
    print('Tier-1 candidate QC complete')
    print('=' * 72)
    print(f'Focal Tier-1 rows: {len(focal)}')
    print(f'Singleton Tier-1 rows: {len(secondary)}')
    print(f'Unique group x CRE candidates: {len(summary)}')
    print(f"Unique candidate CREs: {(summary['dmel_cre_id'].nunique() if not summary.empty else 0)}")
    if not summary.empty:
        print()
        print('Tier-1 scope:')
        print(summary['tier1_scope'].value_counts().to_string())
        print()
        print('Candidate QC:')
        print(summary['candidate_qc'].value_counts().to_string())
        recurrent = summary.loc[summary['candidate_qc'] == 'PASS'].groupby(
            'dmel_cre_id').agg(n_clades=('group_name', 'nunique'), 
                               n_focal_clades=('tier1_scope', lambda x: int((x == 'focal').sum())), 
                               n_secondary_only_clades=('tier1_scope', lambda x: int(
                                   (x == 'secondary_only').sum()))).reset_index().sort_values(['n_clades', 'dmel_cre_id'], ascending=[False, True])
        recurrent = recurrent[recurrent['n_clades'] >= args.recurrence_min_clades]
        if not recurrent.empty:
            print()
            print('QC-passing CREs recurring across clades:')
            print(recurrent.to_string(index=False))
    print()
    print(f'Wrote summary:\n{args.out_summary}')
    print()
    print(f'Wrote details:\n{args.out_detail}')
    print()
    print(f'Wrote metadata:\n{args.metadata_out}')

# ============================================================
# Entry point
# ============================================================

if __name__ == "__main__":
    main()
