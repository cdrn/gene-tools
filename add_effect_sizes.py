#!/usr/bin/env python3
"""
Add effect sizes to all athletic SNPs based on evidence tiers
"""

# Effect size assignments based on evidence strength
EFFECT_SIZES = {
    # Tier 1: Strong evidence (published ORs from meta-analyses)
    'rs1815739': 1.5,   # ACTN3: Meta-analyses show OR ~1.4-1.7 for power
    'rs17602729': 2.17, # AMPD1: OR=2.17 from PMC12152022
    'rs1805086': 2.02,  # MSTN: OR=2.02 from PMC3024427

    # Tier 2: Moderate evidence (well-replicated associations)
    'rs4343': 1.2,      # ACE: Mixed results, conservative estimate
    'rs8192678': 1.3,   # PPARGC1A: Meta-analysis shows consistent association
    'rs4880': 1.3,      # SOD2: 30-40% activity difference
    'rs1800012': 1.3,   # COL1A1: Bone/tendon effects
    'rs35767': 1.3,     # IGF1: T allele higher IGF-1, hypertrophy
    'rs2228570': 1.2,   # VDR FokI: Muscle fiber effects
    'rs1801253': 1.2,   # ADRB1: Cardiac/endurance effects
    'rs2016520': 1.2,   # PPARD: Fat oxidation, endurance enrichment
    'rs6265': 1.2,      # BDNF: Motor learning effects
    'rs1205': 1.2,      # CRP: 20% reduction per allele
    'rs4680': 1.2,      # COMT: Dopamine, stress resilience

    # Tier 3: Preliminary/weaker evidence
    'rs1042713': 1.1,   # ADRB2 Arg16Gly: Inconsistent results
    'rs1042714': 1.0,   # ADRB2 Gln27Glu: Weak performance association
    'rs4994': 1.1,      # ADRB3: Some evidence
    'rs8111989': 1.1,   # CKM: Limited evidence
    'rs2070744': 1.1,   # NOS3: Endurance association
    'rs2010963': 1.1,   # VEGFA: Vascular effects
    'rs11549465': 1.1,  # HIF1A: Rare variant
    'rs1867785': 1.1,   # EPAS1: Ancestry-dependent
    'rs56721780': 1.0,  # EPAS1 Tibetan: Very population-specific
    'rs2228145': 1.1,   # IL6R: Speculative endurance benefit
    'rs1695': 1.1,      # GSTP1: Oxidative stress
    'rs1800795': 1.1,   # IL6: Inflammation
    'rs1800629': 1.1,   # TNF: Inflammation
    'rs12722': 1.1,     # COL5A1: Tendon effects
    'rs1800255': 1.2,   # COL3A1: Injury risk (OR=4.79 for pathology)
    'rs4253778': 1.1,   # PPARA: Fat oxidation
    'rs659366': 1.1,    # UCP2: Energy efficiency
    'rs1049434': 1.1,   # SLC16A1: Lactate transport
    'rs470117': 1.1,    # CPT1B: Beta-oxidation
    'rs9939609': 1.2,   # FTO: OR=1.34-1.55 for obesity (inverse for athleticism)
    'rs6721961': 1.1,   # NFE2L2: Oxidative stress
    'rs1800849': 1.1,   # UCP3: Thermogenesis
    'rs1800497': 1.2,   # DRD2: 40% receptor reduction
}

# Print update SQL-style commands for manual editing
for rsid, effect_size in sorted(EFFECT_SIZES.items()):
    print(f"# {rsid}: effect_size={effect_size}")

print(f"\nTotal SNPs with effect sizes: {len(EFFECT_SIZES)}")
