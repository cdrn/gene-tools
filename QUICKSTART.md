# Quickstart Guide

## 1. Get Your Raw DNA Data

Download your raw DNA file from AncestryDNA:
1. Log into AncestryDNA
2. Go to Settings → Privacy Center → Download Your Data
3. Download the .zip file
4. Extract `AncestryDNA.txt` to this folder

## 2. Install Dependencies

```bash
# Using uv (recommended)
uv pip install snps pandas numpy

# Or using pip
pip install snps pandas numpy
```

## 3. Run Analysis

```bash
# Full health/trait analysis
uv run python analyze_dna.py AncestryDNA.txt

# Athletic performance score only
uv run python athletic_score.py AncestryDNA.txt

# Save output to file
uv run python analyze_dna.py AncestryDNA.txt > my_results.txt
```

## 4. Read Your Results

The analysis will show:
- **Clinical markers** (FDA-recognized drug metabolism variants)
- **Disease risks** (Alzheimer's, diabetes, cardiovascular)
- **Mental health variants** (depression, ADHD, bipolar markers)
- **Physical traits** (eye color, muscle type, caffeine metabolism)
- **Athletic score** (endurance vs power spectrum)

Each variant includes:
- Gene name and rsID
- Your genotype (e.g., AA, AG, GG)
- Interpretation with color coding (red=risk, green=protective)
- Direct link to SNPedia for more info

## Example Output

```
PHARMACOGENOMICS (FDA Recognized):

  CYP2C19 [Clopidogrel (Plavix)] (rs4244285): GG
  Evidence: FDA Black Box Warning
  → Normal metabolizer (*1/*1) - standard effectiveness
  🔗 https://www.snpedia.com/index.php/rs4244285

ATHLETIC PERFORMANCE POLYGENIC SCORE:

Overall Score: +2
Athletic Type: Balanced/Mixed
Performance Spectrum:
Power |●------ Endurance
```

## Other DNA Files

Works with raw data from:
- AncestryDNA (tested)
- 23andMe (should work)
- MyHeritage (should work)
- FamilyTreeDNA (should work)

Just replace `AncestryDNA.txt` with your file name.

## Troubleshooting

**"No module named 'snps'"** → Install dependencies: `uv pip install snps`

**"File not found"** → Make sure your DNA file is in the same folder

**Missing SNPs** → Consumer tests only cover ~0.02% of your genome. Some markers won't be available.