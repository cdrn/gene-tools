# DNA Analysis

Local analysis of AncestryDNA raw data. No data leaves your machine.

## Usage

**This project uses [uv](https://github.com/astral-sh/uv) for Python dependency management.**

1. Place your `AncestryDNA.txt` file in this directory
2. Download GWAS summary statistics (for cognitive score):
   - Register at https://thessgac.com/
   - Download `GWAS_EA_excl23andMe.txt` (Lee et al. 2018)
   - Rename to `EA_GWAS.txt` for easier handling
   - Place in this directory
3. Run the analysis:

```bash
# Full analysis (includes all modules)
uv run analyze_dna.py AncestryDNA.txt

# Run individual modules
uv run athletic_score.py AncestryDNA.txt

# Cognitive score - multiple approaches
uv run cognitive_score.py AncestryDNA.txt                    # Default (p<0.05)
uv run cognitive_score.py AncestryDNA.txt --top 500          # Top 500 SNPs (RECOMMENDED)
uv run cognitive_score.py AncestryDNA.txt --top-analysis     # Compare different top N values
uv run cognitive_score.py AncestryDNA.txt --compare          # Compare p-value thresholds
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

### `cognitive_score.py` - Educational Attainment Polygenic Score
Calculates educational attainment / cognitive ability polygenic score using Lee et al. (2018) GWAS data:
- Based on 766,345 individuals (excluding 23andMe)
- 10.1 million SNPs in full dataset
- Explains ~11-13% of educational attainment variance
- Correlates with intelligence (r ≈ 0.7)
- **NEW: Top SNPs approach** - Use `--top 500` for better accuracy with high-achievers
  - Focuses on highest-impact variants only
  - Reduces noise from weak-effect SNPs
  - Better for detecting exceptional cognitive potential
- Multiple analysis modes:
  - Standard p-value thresholds (0.05, 0.001, 5e-8)
  - Top N SNPs selection (10, 100, 500, 1000, etc.)
  - Comparison modes to find your optimal approach
- **Requires GWAS summary statistics file** (see Usage)

## Data Sources
- FDA Pharmacogenetic Associations
- ACMG Secondary Findings v3.2
- CPIC Guidelines
- Replicated GWAS findings
- SNPedia references for all markers

## Requirements

This project uses **uv** for dependency management. Install uv first:
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Dependencies are managed via `pyproject.toml` and will be automatically installed when you run:
```bash
uv run analyze_dna.py
```

Manual installation (if needed):
```bash
uv pip install snps pandas numpy
```

## Output
- ~300+ SNPs analyzed across clinical and research categories
- Color-coded risk assessment
- Direct SNPedia links for each variant
- Athletic performance score with visual spectrum
- Educational attainment polygenic score (if GWAS data available)
- Y-chromosome haplogroup (for males)