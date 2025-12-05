# DNA Analysis

Local analysis of AncestryDNA raw data. No data leaves your machine.

## Quick Start

**This project uses [uv](https://github.com/astral-sh/uv) for Python dependency management.**

### Simple Command-Line Interface

```bash
# Full analysis (all modules)
uv run python dna_cli.py ancestryDNA.txt

# Athletic Performance Score
uv run python dna_cli.py --athletic ancestryDNA.txt

# Y-chromosome Haplogroup (males only)
uv run python dna_cli.py --haplogroup ancestryDNA.txt

# Help and options
uv run python dna_cli.py --help
```

The CLI includes:
- File validation and helpful error messages
- File size checking to ensure valid DNA files
- Support for AncestryDNA, 23andMe, and other formats
- Automatic format detection

### Advanced Usage

1. Place your `AncestryDNA.txt` file in this directory
2. Run the analysis:

```bash
# Full DNA analysis
uv run python analyze_dna.py AncestryDNA.txt

# Athletic performance score
uv run python athletic_score.py AncestryDNA.txt
```

## Modules

### `analyze_dna.py` - Main Analysis
- **Pharmacogenomics**: Drug metabolism variants (CYP2C9, CYP2C19, CYP2D6, etc.)
- **Disease Risk**: APOE (Alzheimer), Factor V Leiden, hemochromatosis, diabetes (TCF7L2), obesity (FTO)
- **Mental Health**: Serotonin (HTR1A/2A, SLC6A4), dopamine (COMT, DRD2), BDNF, bipolar/schizophrenia markers
- **Physical Traits**: Eye/hair/skin color, muscle type (ACTN3), lactose tolerance, caffeine metabolism
- **Y-Haplogroup**: Paternal lineage prediction

### `athletic_score.py` - Athletic Performance Analysis
Analyzes athletic performance genetics based on validated SNPs:
- Muscle/energy systems (ACTN3, ACE, PPARGC1A, AMPD1)
- Adrenergic response (ADRB2, ADRB3)
- Oxygen/hypoxia (HIF1A, EPAS1)
- Inflammation/recovery (IL6, TNF, COL5A1)
- Fuel metabolism (PPARA, UCP2, SLC16A1)
- Note: Polygenic scoring is under revision for improved accuracy

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
- Athletic performance analysis
- Y-chromosome haplogroup (for males)