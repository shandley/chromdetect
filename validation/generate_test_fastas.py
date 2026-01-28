#!/usr/bin/env python3
"""Generate synthetic FASTA files from NCBI assembly reports for testing ChromDetect."""

import os
import random

def parse_assembly_report(report_path):
    """Parse NCBI assembly report and extract scaffold information."""
    scaffolds = []
    with open(report_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9:
                name = parts[0]
                role = parts[1]
                molecule = parts[2]
                location_type = parts[3]
                genbank = parts[4]
                refseq = parts[6] if len(parts) > 6 else ''
                length = int(parts[8]) if len(parts) > 8 and parts[8].isdigit() else 1000000
                scaffolds.append({
                    'name': name,
                    'role': role,
                    'molecule': molecule,
                    'type': location_type,
                    'genbank': genbank,
                    'refseq': refseq,
                    'length': length
                })
    return scaffolds

def generate_fasta(scaffolds, output_path, use_field='name', max_scaffolds=50):
    """Generate a synthetic FASTA file with realistic-length placeholder sequences."""
    random.seed(42)
    with open(output_path, 'w') as f:
        for i, scaffold in enumerate(scaffolds[:max_scaffolds]):
            name = scaffold.get(use_field, scaffold['name'])
            if not name or name == 'na':
                name = scaffold['name']
            # Use actual length but generate short placeholder sequence
            length = scaffold['length']
            # Write header with length info
            f.write(f">{name}\n")
            # Generate a short representative sequence (100bp per scaffold for testing)
            seq = ''.join(random.choices('ACGT', k=min(100, length)))
            f.write(f"{seq}\n")
    return output_path

def main():
    reports = [
        # Original assemblies
        ('zebra_finch_report.txt', 'zebra_finch', 'SUPER_* naming (VGP)'),
        ('chicken_report.txt', 'chicken', 'Numeric naming (GRC)'),
        ('human_grch38_report.txt', 'human', 'Numeric naming (GRCh38)'),
        ('zebrafish_report.txt', 'zebrafish', 'Numeric naming (GRCz11)'),
        # Expanded validation assemblies
        ('arabidopsis_report.txt', 'arabidopsis', 'Numeric naming (TAIR10 plant)'),
        ('ecoli_report.txt', 'ecoli', 'Accession naming (bacterial)'),
        ('t2t_human_report.txt', 't2t_human', 'Numeric naming (T2T-CHM13)'),
    ]

    for report_file, prefix, description in reports:
        if not os.path.exists(report_file):
            print(f"Skipping {report_file} - not found")
            continue

        print(f"\n{'='*60}")
        print(f"Processing: {prefix} ({description})")
        print(f"{'='*60}")

        scaffolds = parse_assembly_report(report_file)
        print(f"Found {len(scaffolds)} scaffolds in assembly report")

        # Count by role
        roles = {}
        for s in scaffolds:
            role = s['role']
            roles[role] = roles.get(role, 0) + 1
        print(f"Roles: {roles}")

        # Generate FASTA using scaffold names
        fasta_path = f"{prefix}_test.fasta"
        generate_fasta(scaffolds, fasta_path, use_field='name')
        print(f"Generated: {fasta_path}")

        # Also generate using RefSeq accessions for variety
        refseq_path = f"{prefix}_refseq.fasta"
        generate_fasta(scaffolds, refseq_path, use_field='refseq')
        print(f"Generated: {refseq_path}")

if __name__ == '__main__':
    main()
