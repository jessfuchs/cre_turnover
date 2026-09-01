#!/usr/bin/env python3

# ============================================================
# Assign candidate genes to Tier-1 CREs
#
# Purpose:
#   Assign primary and secondary candidate genes to prioritized
#   Tier-1 CREs using D. melanogaster flanking-gene annotations.
#
# Assignment:
#   Candidate genes are ranked by absolute CRE-to-gene
#   distance. Equal-distance ties are resolved deterministically.
#   Agreement with the original SCRMshaw target-gene annotation
#   is retained as QC information.
#
# Input/output paths:
#   Supplied by the pipeline wrapper using
#   config/candidate_config.sh.
# ============================================================

import argparse
from datetime import datetime
import hashlib
from pathlib import Path
import platform
import sys

import pandas as pd

# ============================================================
# Arguments
# ============================================================

def parse_args():
    parser = argparse.ArgumentParser(description='Assign primary and secondary candidate genes to Tier-1 CRE candidates using D. melanogaster flanking-gene distances.')
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--metadata-out', type=Path, required=True)
    return parser.parse_args()

# ============================================================
# Required columns
# ============================================================

REQUIRED_COLUMNS = {
    "dmel_cre_id",
    "candidate_priority",
    "downstream_priority",
    "fbgn_target_genes",

    "dmel_fbgn_flanking_gene",
    "distance_flanking_gene",

    "dmel_fbgn_next_flanking_gene",
    "distance_next_gene",

    "gene_distance_qc",
}

# ============================================================
# Helpers
# ============================================================

def require_file(path):
    if not path.is_file():
        raise SystemExit(f'ERROR: required input file not found:\n{path}')

def require_columns(df, required, label):
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise SystemExit(f'ERROR: required columns missing from {label}:\n' + '\n'.join(missing))

def normalize_missing(value):
    if pd.isna(value):
        return 'NA'
    value = str(value).strip()
    if value.lower() in {'', 'na', 'nan', 'none'}:
        return 'NA'
    return value

def parse_target_genes(value):
    """
    Convert pipe-separated FBgn target genes into a set.
    """
    return set(parse_fbgns(value))

def parse_fbgns(value):
    """
    Parse pipe-separated FBgn identifiers.
    """
    raw = normalize_missing(value)
    if raw == 'NA':
        return []
    return sorted({gene.strip() for gene in raw.split('|') if gene.strip() and gene.strip().lower() not in {'na', 'nan', 'none'}})

def file_sha256(path):
    sha = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            sha.update(block)
    return sha.hexdigest()

def classify_relation(distance):
    """
    Describe CRE-gene relationship based on distance.

    distance == 0
        overlapping

    distance > 0
        nearest_nonoverlapping

    missing
        unknown
    """
    if pd.isna(distance):
        return 'unknown'
    if float(distance) == 0:
        return 'overlapping'
    return 'nearest_nonoverlapping'

# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()
    require_file(args.input)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.metadata_out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.input, sep='\t', dtype=str).fillna('')
    require_columns(df, REQUIRED_COLUMNS, 'candidate gene-distance table')

    # --------------------------------------------------------
    # One row per CRE
    # --------------------------------------------------------

    if df['dmel_cre_id'].duplicated().any():
        duplicated = df.loc[df['dmel_cre_id'].duplicated(keep=False), 'dmel_cre_id'].drop_duplicates().tolist()
        raise SystemExit('ERROR: duplicate CRE IDs in input:\n' + '\n'.join(sorted(duplicated)))

    # ========================================================
    # Numeric distances
    # ========================================================

    for column in ['distance_flanking_gene', 'distance_next_gene']:
        df[column] = pd.to_numeric(df[column], errors='coerce')

    # ========================================================
    # Assign genes
    # ========================================================

    assignment_rows = []
    for _, row in df.iterrows():
        cre_id = row['dmel_cre_id']
        target_genes = parse_target_genes(row['fbgn_target_genes'])

        # ----------------------------------------------------
        # Raw distance-derived gene annotations
        # ----------------------------------------------------
        
        fbgn1 = normalize_missing(row['dmel_fbgn_flanking_gene'])
        fbgn2 = normalize_missing(row['dmel_fbgn_next_flanking_gene'])
        fbgn1_values = parse_fbgns(fbgn1)
        fbgn2_values = parse_fbgns(fbgn2)
        d1 = row['distance_flanking_gene']
        d2 = row['distance_next_gene']
        
        # ====================================================
        # Build distance-derived candidate genes
        #
        # Distances are ranked by absolute genomic distance.
        # Raw signed distances are retained separately in the
        # provenance columns of the output.
        # ====================================================
        
        candidates = []
        if fbgn1_values and (not pd.isna(d1)):
            for fbgn in fbgn1_values:
                candidates.append({'fbgn': fbgn, 'distance': abs(float(d1)), 'source': 'flanking_gene'})
        if fbgn2_values and (not pd.isna(d2)):
            for fbgn in fbgn2_values:
                candidates.append({'fbgn': fbgn, 'distance': abs(float(d2)), 'source': 'next_flanking_gene'})

        # ----------------------------------------------------
        # Deduplicate identical FBgn assignments
        #
        # If the same FBgn appears twice, retain its shortest
        # distance and one deterministic source label.
        # ----------------------------------------------------

        collapsed = {}
        for candidate in candidates:
            fbgn = candidate['fbgn']
            if fbgn not in collapsed:
                collapsed[fbgn] = candidate
            else:
                existing = collapsed[fbgn]
                if candidate['distance'] < existing['distance']:
                    collapsed[fbgn] = candidate
        candidates = list(collapsed.values())


        # ----------------------------------------------------
        # Stable deterministic ranking
        # ----------------------------------------------------

        candidates = sorted(candidates, key=lambda x: (x['distance'], x['fbgn'], x['source']))

        # ====================================================
        # Primary / secondary assignment
        # ====================================================

        if len(candidates) == 0:
            primary_fbgn = 'NA'
            primary_distance = pd.NA
            primary_source = 'NA'
            secondary_fbgn = 'NA'
            secondary_distance = pd.NA
            secondary_source = 'NA'
        elif len(candidates) == 1:
            primary_fbgn = candidates[0]['fbgn']
            primary_distance = candidates[0]['distance']
            primary_source = candidates[0]['source']
            secondary_fbgn = 'NA'
            secondary_distance = pd.NA
            secondary_source = 'NA'
        else:
            primary_fbgn = candidates[0]['fbgn']
            primary_distance = candidates[0]['distance']
            primary_source = candidates[0]['source']
            secondary_fbgn = candidates[1]['fbgn']
            secondary_distance = candidates[1]['distance']
            secondary_source = candidates[1]['source']


        # ====================================================
        # Relations
        # ====================================================

        primary_relation = classify_relation(primary_distance)
        secondary_relation = classify_relation(secondary_distance)

        # ====================================================
        # Gene assignment QC
        # ====================================================

        qc_flags = []

        distance_qc = normalize_missing(row['gene_distance_qc'])
        if distance_qc != 'PASS':
            qc_flags.append(f'GENE_DISTANCE_QC_{distance_qc}')
    
        if primary_fbgn == 'NA':
            qc_flags.append('NO_DISTANCE_DERIVED_GENE')

        # ----------------------------------------------------
        # Check against original SCRMshaw FBgn targets
        # ----------------------------------------------------

        if primary_fbgn != 'NA' and primary_fbgn not in target_genes:
            qc_flags.append('PRIMARY_FBGN_NOT_IN_REFERENCE_TARGETS')
            
        if secondary_fbgn != 'NA' and secondary_fbgn not in target_genes:
            qc_flags.append('SECONDARY_FBGN_NOT_IN_REFERENCE_TARGETS')
            
        if not target_genes:
            qc_flags.append('NO_REFERENCE_TARGET_GENES')

        # ----------------------------------------------------
        # Equal distance tie
        # ----------------------------------------------------

        distance_tie = primary_fbgn != 'NA' and secondary_fbgn != 'NA' and (not pd.isna(primary_distance)) and (not pd.isna(secondary_distance)) and (float(primary_distance) == float(secondary_distance))
        if distance_tie:
            qc_flags.append('PRIMARY_SECONDARY_DISTANCE_TIE')

        # ----------------------------------------------------
        # Determine whether primary gene is one of original
        # SCRMshaw target-gene assignments
        # ----------------------------------------------------

        primary_in_reference = primary_fbgn != 'NA' and primary_fbgn in target_genes
        secondary_in_reference = secondary_fbgn != 'NA' and secondary_fbgn in target_genes

        # ----------------------------------------------------
        # Final QC
        # ----------------------------------------------------

        qc_flags = sorted(set(qc_flags))
        assignment_qc = 'PASS' if not qc_flags else ';'.join(qc_flags)

        # ====================================================
        # Output row
        # ====================================================

        assignment_rows.append({'dmel_cre_id': cre_id, 
                                'candidate_priority': row['candidate_priority'], 
                                'downstream_priority': row['downstream_priority'], 
                                'n_total_tier1_clades': normalize_missing(row.get('n_total_tier1_clades', 'NA')), 
                                'n_focal_tier1_clades': normalize_missing(row.get('n_focal_tier1_clades', 'NA')), 
                                'n_secondary_only_clades': normalize_missing(row.get('n_secondary_only_clades', 'NA')), 
                                'tier1_clades': normalize_missing(row.get('tier1_clades', 'NA')), 
                                'fbgn_target_genes': normalize_missing(row['fbgn_target_genes']), 
                                'n_reference_target_genes': len(target_genes), 
                                'primary_candidate_fbgn': primary_fbgn, 
                                'primary_candidate_distance_bp': primary_distance, 
                                'primary_candidate_relation': primary_relation, 
                                'primary_candidate_source': primary_source, 
                                'primary_in_reference_targets': 'yes' if primary_in_reference else 'no', 
                                'secondary_candidate_fbgn': secondary_fbgn, 
                                'secondary_candidate_distance_bp': secondary_distance, 
                                'secondary_candidate_relation': secondary_relation, 
                                'secondary_candidate_source': secondary_source, 
                                'secondary_in_reference_targets': 'yes' if secondary_in_reference else 'no', 
                                'primary_secondary_distance_tie': 'yes' if distance_tie else 'no', 
                                'dmel_fbgn_flanking_gene': fbgn1, 
                                'distance_flanking_gene': d1, 
                                'dmel_fbgn_next_flanking_gene': fbgn2, 
                                'distance_next_gene': d2, 
                                'gene_distance_qc': distance_qc, 
                                'gene_assignment_qc': assignment_qc})

    # ========================================================
    # Build output
    # ========================================================

    out = pd.DataFrame(assignment_rows)

    # --------------------------------------------------------
    # Stable ordering
    # --------------------------------------------------------

    priority_order = {'high': 1, 
                      'medium': 2, 
                      'exploratory': 3}

    out['_priority_rank'] = out['downstream_priority'].map(priority_order).fillna(99)

    out = out.sort_values(['_priority_rank', 'dmel_cre_id'], kind='mergesort').drop(columns=['_priority_rank']).reset_index(drop=True)

    # --------------------------------------------------------
    # Final uniqueness QC
    # --------------------------------------------------------

    if out['dmel_cre_id'].duplicated().any():
        raise SystemExit('ERROR: duplicate CRE IDs in gene-assignment output.')

    # ========================================================
    # Write
    # ========================================================

    out.to_csv(args.out, sep="\t", index=False)

    # ========================================================
    # Metadata
    # ========================================================

    metadata = pd.DataFrame([{'script': Path(__file__).name, 
                              'run_timestamp': datetime.now().astimezone().isoformat(), 
                              'python_version': sys.version.split()[0], 
                              'pandas_version': pd.__version__, 
                              'platform': platform.platform(), 
                              'input_file': str(args.input.resolve()), 
                              'input_sha256': file_sha256(args.input), 
                              'n_candidate_cres': len(out), 
                              'n_primary_assigned': int((out['primary_candidate_fbgn'] != 'NA').sum()), 
                              'n_secondary_assigned': int((out['secondary_candidate_fbgn'] != 'NA').sum()), 
                              'n_assignment_pass': int((out['gene_assignment_qc'] == 'PASS').sum()), 
                              'n_primary_in_reference_targets': int((out['primary_in_reference_targets'] == 'yes').sum()), 
                              'n_secondary_in_reference_targets': int((out['secondary_in_reference_targets'] == 'yes').sum()), 
                              'n_distance_ties': int((out['primary_secondary_distance_tie'] == 'yes').sum())}])
    
    metadata.to_csv(args.metadata_out, sep='\t', index=False)

    # ========================================================
    # Console summary
    # ========================================================

    print()
    print('=' * 72)
    print('Candidate gene assignment complete')
    print('=' * 72)
    print(f'Candidate CREs: {len(out)}')
    print(f"Primary candidate genes assigned: {(out['primary_candidate_fbgn'] != 'NA').sum()}")
    print(f"Secondary candidate genes assigned: {(out['secondary_candidate_fbgn'] != 'NA').sum()}")
    print()
    print('Gene assignment QC:')
    print(out['gene_assignment_qc'].value_counts().to_string())
    flagged = out.loc[out['gene_assignment_qc'] != 'PASS']
    if not flagged.empty:
        print()
        print('Assignments requiring inspection:')
        print(flagged[['dmel_cre_id', 
                       'candidate_priority', 
                       'fbgn_target_genes', 
                       'primary_candidate_fbgn', 
                       'primary_candidate_distance_bp', 
                       'primary_in_reference_targets', 
                       'secondary_candidate_fbgn', 
                       'secondary_candidate_distance_bp', 
                       'secondary_in_reference_targets', 
                       'gene_assignment_qc']].to_string(index=False))
    print()
    print(f'Wrote gene assignments:\n{args.out}')
    print()
    print(f'Wrote metadata:\n{args.metadata_out}')

# ============================================================
# Entry point
# ============================================================

if __name__ == "__main__":
    main()
