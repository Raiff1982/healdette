Results Note
Date: September 12, 2025
Version: 1.0.0
Commit: main-2025-09-12
DOI: 10.57967/hf/5917

Execution Details:
```bash
python run_pipeline.py --deterministic --seed 42
```

Environment:
- Python 3.10.8
- Environment hash: <SHA256 of pip freeze output>
- OS: Windows 10
- Hardware: CPU-only execution

Input Parameters:
- Ancestry profile: Native, Irish
- HLA matches: 2
- Prior exposures: SARS-CoV-2, Influenza-B
- Metabolic factor: 1.2
- Random seed: 42

Generated Sequences Analysis:
| ID | Length | Score | Disorder | Cys Pairs | N-glyc | GRAVY | pI |
|----|--------|--------|----------|------------|--------|--------|-----|

Key Findings:
1. Length distribution: 43-563 amino acids
2. Personalization scores: 0.62-0.78
3. Disorder scores: 0.185-0.677
4. Glycosylation sites: 7 total (avg 1.0 per sequence)
5. Cysteine pairs: 3/7 sequences have paired cysteines

Validation Status:
- Environment: See environment.yaml
- Checksums: See checksums.sha256
- Full results: validation_results_20250912_152239.json

For reproduction:
1. Clone repository
2. Install dependencies from environment.yaml
3. Run: python run_pipeline.py --deterministic
4. Verify checksums