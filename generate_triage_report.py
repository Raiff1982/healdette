"""
Generate detailed triage report for antibody designs.
"""

import pandas as pd
import json
from datetime import datetime

def create_triage_report(results_json, output_file):
    """Create a detailed triage report in markdown format."""
    with open(results_json, 'r') as f:
        data = json.load(f)
    
    report = []
    report.append("# Antibody Design Triage Report")
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    # Summary statistics
    report.append("## Summary Statistics")
    report.append("| Metric | Value |")
    report.append("| --- | --- |")
    report.append(f"| Total Sequences | {len(data['personalized_binders'])} |")
    
    # Triage table
    report.append("\n## Sequence Analysis")
    report.append("| ID | Length | Score | Disorder | Cys Pairs | Glyco | GRAVY | Status |")
    report.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
    
    for i, binder in enumerate(data['validated_binders'], 1):
        val = binder['validation']
        status = "PASS" if (
            val['disorder'] <= 0.5 and
            not val['signal_peptide']['has_signal'] and
            val['cysteines']['paired'] and
            -1.0 <= val['properties']['GRAVY'] <= 1.0
        ) else "FAIL"
        
        report.append(
            f"| {i} | {len(binder['sequence'])} | "
            f"{binder['personalization_score']:.3f} | "
            f"{val['disorder']:.3f} | "
            f"{val['cysteines']['count']//2} | "
            f"{len(val['glycosylation'])} | "
            f"{val['properties']['GRAVY']:.3f} | "
            f"{status} |"
        )
    
    # Failure analysis
    report.append("\n## Failure Analysis")
    failure_counts = {
        "High Disorder": sum(1 for b in data['validated_binders'] 
                           if b['validation']['disorder'] > 0.5),
        "Signal Peptide": sum(1 for b in data['validated_binders'] 
                            if b['validation']['signal_peptide']['has_signal']),
        "Unpaired Cys": sum(1 for b in data['validated_binders'] 
                           if not b['validation']['cysteines']['paired']),
        "GRAVY Outside Range": sum(1 for b in data['validated_binders'] 
                                 if not -1.0 <= b['validation']['properties']['GRAVY'] <= 1.0)
    }
    
    for reason, count in failure_counts.items():
        report.append(f"- {reason}: {count} sequences")
    
    # Write report
    with open(output_file, 'w') as f:
        f.write('\n'.join(report))

if __name__ == "__main__":
    results_json = "output/validation_results_20250912_152239.json"
    output_file = "output/triage_report.md"
    create_triage_report(results_json, output_file)