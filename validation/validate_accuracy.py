#!/usr/bin/env python3
"""
Validate ChromDetect classification accuracy against NCBI assembly reports.

This script compares ChromDetect's scaffold classifications against the
ground truth from NCBI assembly reports.
"""

import json
import os
import sys

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from chromdetect import classify_fasta


def parse_assembly_report(report_path):
    """Parse NCBI assembly report and extract ground truth classifications."""
    ground_truth = {}
    with open(report_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9:
                name = parts[0]
                role = parts[1]
                refseq = parts[6] if len(parts) > 6 else ''

                # Map NCBI roles to ChromDetect classifications
                if role == 'assembled-molecule':
                    classification = 'chromosome'
                elif role == 'unlocalized-scaffold':
                    classification = 'unlocalized'
                else:  # unplaced-scaffold, alt-scaffold, fix-patch, novel-patch
                    classification = 'unplaced'

                ground_truth[name] = classification
                if refseq and refseq != 'na':
                    ground_truth[refseq] = classification

    return ground_truth


def validate_assembly(fasta_path, report_path, description):
    """Validate ChromDetect against ground truth from assembly report."""
    print(f"\n{'='*70}")
    print(f"Validation: {description}")
    print(f"FASTA: {fasta_path}")
    print(f"Report: {report_path}")
    print(f"{'='*70}")

    # Parse ground truth
    ground_truth = parse_assembly_report(report_path)

    # Run ChromDetect
    results, stats = classify_fasta(fasta_path)

    # Compare classifications
    correct = 0
    incorrect = 0
    not_found = 0

    details = {'correct': [], 'incorrect': [], 'not_found': []}

    for scaffold in results:
        name = scaffold.name
        predicted = scaffold.classification

        if name in ground_truth:
            expected = ground_truth[name]
            if predicted == expected:
                correct += 1
                details['correct'].append((name, predicted))
            else:
                incorrect += 1
                details['incorrect'].append((name, predicted, expected))
        else:
            not_found += 1
            details['not_found'].append((name, predicted))

    total = correct + incorrect
    accuracy = (correct / total * 100) if total > 0 else 0

    print(f"\nResults:")
    print(f"  Total scaffolds classified: {len(results)}")
    print(f"  Matched to ground truth: {total}")
    print(f"  Correct: {correct}")
    print(f"  Incorrect: {incorrect}")
    print(f"  Not in report: {not_found}")
    print(f"  Accuracy: {accuracy:.1f}%")

    if details['incorrect']:
        print(f"\nMisclassifications (first 10):")
        for name, predicted, expected in details['incorrect'][:10]:
            print(f"  {name}: predicted={predicted}, expected={expected}")

    return {
        'description': description,
        'fasta': fasta_path,
        'total': len(results),
        'matched': total,
        'correct': correct,
        'incorrect': incorrect,
        'accuracy': accuracy,
        'misclassifications': details['incorrect'][:10]
    }


def main():
    """Run validation on all test assemblies."""
    validations = [
        ('zebra_finch_test.fasta', 'zebra_finch_report.txt',
         'Zebra Finch VGP (SUPER_* naming)'),
        ('zebra_finch_refseq.fasta', 'zebra_finch_report.txt',
         'Zebra Finch VGP (RefSeq NC_* naming)'),
        ('chicken_test.fasta', 'chicken_report.txt',
         'Chicken GRC (Numeric naming)'),
        ('human_test.fasta', 'human_grch38_report.txt',
         'Human GRCh38 (Numeric naming)'),
        ('zebrafish_test.fasta', 'zebrafish_report.txt',
         'Zebrafish GRCz11 (Numeric naming)'),
    ]

    results = []
    for fasta, report, description in validations:
        if os.path.exists(fasta) and os.path.exists(report):
            result = validate_assembly(fasta, report, description)
            results.append(result)
        else:
            print(f"Skipping {description} - files not found")

    # Summary
    print(f"\n{'='*70}")
    print("VALIDATION SUMMARY")
    print(f"{'='*70}")
    print(f"\n{'Assembly':<45} {'Accuracy':<10} {'Correct/Total':<15}")
    print("-" * 70)

    total_correct = 0
    total_matched = 0

    for r in results:
        desc = r['description'][:44]
        print(f"{desc:<45} {r['accuracy']:>6.1f}%    {r['correct']}/{r['matched']}")
        total_correct += r['correct']
        total_matched += r['matched']

    overall_accuracy = (total_correct / total_matched * 100) if total_matched > 0 else 0
    print("-" * 70)
    print(f"{'OVERALL':<45} {overall_accuracy:>6.1f}%    {total_correct}/{total_matched}")

    # Save results to JSON
    with open('validation_results.json', 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nDetailed results saved to: validation_results.json")


if __name__ == '__main__':
    main()
