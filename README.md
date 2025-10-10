# DNA Analysis Script

Comprehensive analysis of your AncestryDNA genetic data **100% offline** - no data transmitted to remote servers.

## Usage

```bash
source .venv/bin/activate
python analyze_dna.py
```

## What It Analyzes

### Tier 1: Clinical Grade Markers
- **Pharmacogenomics**: Drug metabolism (warfarin, clopidogrel, SSRIs, codeine, etc.)
- **Cardiovascular**: APOE (Alzheimer risk), Factor V Leiden, hemochromatosis
- **Blood Clotting**: Prothrombin, clotting disorders
- Sources: FDA Table, ACMG Secondary Findings, CPIC Guidelines

### Health & Disease Risk
- **Metabolic**: Type 2 diabetes (TCF7L2, PPARG), obesity (FTO), fatty liver
- **Cardiovascular**: CAD risk, triglycerides, inflammation markers
- **Cancer**: Tumor suppressor genes, immune markers
- **Autoimmune**: RA, T1D, celiac disease risk
- **Bone & Kidney**: Osteoporosis, kidney disease markers

### Mood & Mental Health (Extensive)
- **Serotonin**: HTR1A, HTR2A, TPH2, SLC6A4 (depression, anxiety, SSRI response)
- **Dopamine**: COMT (Warrior/Worrier), DRD2, DRD4, DBH (ADHD, reward, addiction)
- **Neuroplasticity**: BDNF (learning, memory, depression)
- **Bipolar & Schizophrenia**: CACNA1C, ANK3, ZNF804A, NRG1, TCF4
- **Anxiety & Stress**: FAAH ("happiness gene"), OXTR (empathy, social bonding)
- **Other**: Glutamate receptors, opioid receptors, folate/mood connection

### Physical Traits
- **Appearance**: Eye color, hair color, skin color, freckles, hair texture
- **Body**: Muscle fiber type (ACTN3), height, earwax type, body odor
- **Metabolism**: Lactose tolerance, bitter taste, alcohol flush, caffeine metabolism
- **Substances**: Nicotine dependence, alcohol metabolism
- **Senses**: Cilantro taste, asparagus odor detection
- **Sleep**: Circadian rhythm, chronotype, sleep quality
- **Vitamins**: Vitamin D, B12 absorption

### Ancestry
- **Y-Haplogroup**: Paternal lineage (males only) - informational, limited accuracy
- **mtDNA**: Maternal lineage overview

## Features

- **200+ SNP markers** across all categories
- **SNPedia links** for every marker
- **Evidence citations** for clinical markers
- **Clean, readable output** without excessive warnings
- **100% local processing** - your DNA never leaves your computer

## Your Results Summary

Based on earlier analysis of your AncestryDNA data:

### Basic Stats
- **Total SNPs**: 677,428
- **Data Completeness**: 99.46%
- **Genome Build**: GRCh37 (Build 37)
- **Sex**: Male

### Y-Haplogroup (Paternal Line)
- **Predicted Haplogroup**: I
- Confidence: 40% (based on 2 markers)
- European origin, likely Scandinavian or Southeastern European
- For more detailed analysis: YFull.com, FamilyTreeDNA, YSEQ.net

### Notable Markers (from previous analysis)
- **Eye Color**: Brown eyes (rs12913832: AA)
- **Hair Color**: Light hair tendency (rs12896399: TT)
- **Lactose Tolerance**: Tolerant (rs4988235: AA)
- **Muscle Type**: Endurance-oriented (rs1815739: TT - ACTN3 "non-sprinter")
- **Caffeine**: Fast metabolizer (rs762551: AA)
- **COMT**: Met/Met "Worrier" type (rs4680: AA)
- **Oxytocin**: High empathy (rs53576: GG)
- **Cardiovascular**: Elevated CAD risk marker (rs1333049)

## Privacy & Security

✓ All processing is done **100% locally** on your machine
✓ **No network calls** are made by this script
✓ **No data** is sent to remote servers
✓ Your genetic data stays on your computer

## Requirements

- Python 3.x
- snps package (already installed via uv)
- pandas, numpy (dependencies of snps)

## Evidence Sources

- **FDA Table of Pharmacogenetic Associations**
- **ACMG Secondary Findings v3.2**
- **CPIC Guidelines** (Clinical Pharmacogenetics Implementation Consortium)
- **GWAS Catalog** (genome-wide association studies)
- **Published research studies** (replicated findings)

## Further Analysis

For more specialized analysis:

- **Y-DNA Haplogroup**: YFull, FamilyTreeDNA, YSEQ
- **mtDNA Haplogroup**: James Lick mtDNA, FamilyTreeDNA
- **Ethnicity/Ancestry**: GEDmatch, DNA.land
- **Genealogy**: Ancestry.com, 23andMe, MyHeritage

## Notes

- Consumer DNA tests (like AncestryDNA) analyze ~677,000 SNPs, not your whole genome
- Genetic risk is probabilistic, not deterministic
- Many traits are influenced by multiple genes and environmental factors
- Haplogroup predictions from consumer tests are rough estimates
- For clinical decisions, consult healthcare professionals

## License

This script is provided as-is for personal use. Use at your own risk.
