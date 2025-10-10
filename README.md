# DNA Analysis

Local analysis of AncestryDNA raw data. No data leaves your machine.

## Usage

1. Place your `AncestryDNA.txt` file in this directory
2. Run the analysis:

```bash
uv run python analyze_dna.py AncestryDNA.txt

# Run athletic performance score separately
uv run python athletic_score.py AncestryDNA.txt
```

## Modules

### `analyze_dna.py` - Main Analysis
- **Pharmacogenomics**: Drug metabolism variants (CYP2C9, CYP2C19, CYP2D6, etc.)
- **Disease Risk**: APOE (Alzheimer), Factor V Leiden, hemochromatosis, diabetes (TCF7L2), obesity (FTO)
- **Mental Health**: Serotonin (HTR1A/2A, SLC6A4), dopamine (COMT, DRD2), BDNF, bipolar/schizophrenia markers
- **Physical Traits**: Eye/hair/skin color, muscle type (ACTN3), lactose tolerance, caffeine metabolism
- **Y-Haplogroup**: Paternal lineage prediction

### `athletic_score.py` - Athletic Polygenic Score
Calculates endurance vs power spectrum based on 19 validated SNPs:
- Muscle/energy systems (ACTN3, ACE, PPARGC1A, AMPD1)
- Adrenergic response (ADRB2, ADRB3)
- Oxygen/hypoxia (HIF1A, EPAS1)
- Inflammation/recovery (IL6, TNF, COL5A1)
- Fuel metabolism (PPARA, UCP2, SLC16A1)

## Data Sources
- FDA Pharmacogenetic Associations
- ACMG Secondary Findings v3.2
- CPIC Guidelines
- Replicated GWAS findings
- SNPedia references for all markers

## Requirements
```bash
pip install snps pandas numpy
# or
uv pip install snps pandas numpy
```

## Output
- ~300+ SNPs analyzed across clinical and research categories
- Color-coded risk assessment
- Direct SNPedia links for each variant
- Athletic performance score with visual spectrum