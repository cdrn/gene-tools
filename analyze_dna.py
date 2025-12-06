#!/usr/bin/env python3
"""
Comprehensive DNA Analysis - All-in-One Script
Analyzes SNP data for health, traits, ancestry, mood, and more
100% local processing - no data transmitted
"""

import snps
import sys
from collections import defaultdict
from athletic_score import analyze_athletic_performance


# ANSI color codes
class Colors:
    HEADER = '\033[95m'      # Magenta
    BLUE = '\033[94m'        # Blue
    CYAN = '\033[96m'        # Cyan
    GREEN = '\033[92m'       # Green
    YELLOW = '\033[93m'      # Yellow
    RED = '\033[91m'         # Red
    BOLD = '\033[1m'         # Bold
    UNDERLINE = '\033[4m'    # Underline
    END = '\033[0m'          # Reset
    GRAY = '\033[90m'        # Gray
    ORANGE = '\033[38;5;208m'  # Orange


def print_section(title):
    """Print a formatted section header"""
    print(f"\n{Colors.BOLD}{Colors.CYAN}{'='*80}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN} {title}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}{'='*80}{Colors.END}\n")


# ============================================================================
# TIER 1: CLINICAL GRADE MARKERS
# FDA recognized, ACMG reportable, strong clinical evidence
# ============================================================================

TIER1_PHARMACOGENOMICS = {
    'rs1799853': {  # CYP2C9*2
        'gene': 'CYP2C9', 'drug': 'Warfarin',
        'evidence': 'FDA Table, CPIC Level A',
        'CC': 'Normal metabolism (*1/*1) - standard dosing',
        'CT': 'Intermediate (*1/*2) - ~20% reduced metabolism, may need 15-30% lower dose',
        'TT': 'Reduced metabolism (*2/*2) - ~40% reduction, may need 40-50% lower dose',
    },
    'rs1057910': {  # CYP2C9*3
        'gene': 'CYP2C9', 'drug': 'Warfarin',
        'evidence': 'FDA Table, CPIC Level A',
        'AA': 'Normal metabolism (*1/*1)',
        'AC': 'Intermediate (*1/*3) - lower dose needed',
        'CC': 'Poor metabolism (*3/*3) - much lower dose',
    },
    'rs4244285': {  # CYP2C19*2
        'gene': 'CYP2C19', 'drug': 'Clopidogrel (Plavix)',
        'evidence': 'FDA Black Box Warning',
        'GG': 'Normal metabolizer (*1/*1) - standard effectiveness',
        'AG': 'Intermediate (*1/*2) - reduced effectiveness',
        'AA': 'Poor metabolizer (*2/*2) - may need alternative drug',
    },
    'rs4986893': {  # CYP2C19*3
        'gene': 'CYP2C19', 'drug': 'Clopidogrel',
        'evidence': 'FDA Table, CPIC',
        'GG': 'Normal (*1/*1)',
        'GA': 'Intermediate (*1/*3)',
        'AA': 'Poor metabolizer (*3/*3)',
    },
    'rs12248560': {  # CYP2C19*17
        'gene': 'CYP2C19', 'drug': 'SSRIs, PPIs',
        'evidence': 'CPIC - ultra-rapid metabolism',
        'CC': 'Normal metabolism',
        'CT': 'Increased activity (*1/*17)',
        'TT': 'Ultra-rapid (*17/*17) - higher SSRI side effects',
    },
    'rs1065852': {  # CYP2D6*4
        'gene': 'CYP2D6', 'drug': 'Codeine, antidepressants, antipsychotics',
        'evidence': 'FDA Table - 25% of drugs',
        'CC': 'Normal metabolizer',
        'CT': 'Intermediate',
        'TT': 'Poor metabolizer (*4/*4) - high drug sensitivity',
    },
    'rs3745274': {  # CYP2B6*6
        'gene': 'CYP2B6', 'drug': 'Efavirenz (HIV)',
        'evidence': 'FDA Table, CPIC',
        'GG': 'Normal metabolizer',
        'GT': 'Intermediate (*1/*6) - may have altered efavirenz response',
        'TT': 'Slow metabolizer (*6/*6) - reduced enzyme activity, higher side effects, lower doses needed',
    },
    'rs1801133': {  # MTHFR C677T
        'gene': 'MTHFR', 'drug': 'Methotrexate, 5-FU chemo',
        'evidence': 'FDA Table, CPIC',
        'CC': 'Normal enzyme function',
        'CT': 'Reduced function (~65% activity)',
        'TT': 'Low function (~10-20%) - need higher folate, dose adjustment',
    },
    'rs1801131': {  # MTHFR A1298C
        'gene': 'MTHFR', 'drug': 'Folate metabolism',
        'evidence': 'FDA Table, CPIC',
        'AA': 'Normal function',
        'AC': 'Slightly reduced (~80% activity)',
        'CC': 'Mildly reduced function (~60% activity)',
    },
    'rs1045642': {  # ABCB1 C3435T
        'gene': 'ABCB1', 'drug': 'Many drugs (P-glycoprotein)',
        'evidence': 'Affects drug transport',
        'CC': 'Normal drug metabolism, cannabis dependence risk, lower cancer risk',
        'CT': 'Slower metaboliser for some drugs',
        'TT': 'Altered drug metabolism/bioavailability, moderately increased cancer risk',
    },
    'rs1142345': {  # TPMT*3C
        'gene': 'TPMT', 'drug': 'Azathioprine, 6-mercaptopurine (immunosuppressants)',
        'evidence': 'FDA Table, CPIC Level A',
        'AA': 'Normal TPMT activity',
        'AG': 'Intermediate activity - lower dose needed',
        'GG': 'Low activity - much lower dose or alternative',
    },
    'rs1801280': {  # NAT2*5
        'gene': 'NAT2', 'drug': 'Isoniazid, hydralazine, sulfonamides',
        'evidence': 'FDA Table - affects acetylation speed',
        'TT': 'Normal/rapid acetylator',
        'TC': 'Intermediate acetylator',
        'CC': 'Slow acetylator (higher drug side effects)',
    },
    'rs4149056': {  # SLCO1B1
        'gene': 'SLCO1B1', 'drug': 'Statins (simvastatin)',
        'evidence': 'FDA Table, CPIC - myopathy risk',
        'TT': 'Normal function',
        'TC': 'Intermediate function - higher myopathy risk',
        'CC': 'Reduced function - high myopathy risk (~17x)',
    },
    'rs9923231': {  # VKORC1
        'gene': 'VKORC1', 'drug': 'Warfarin',
        'evidence': 'FDA Table, CPIC - primary warfarin sensitivity',
        'CC': 'Normal warfarin sensitivity - standard dosing',
        'CT': 'Increased sensitivity - reduced warfarin dose needed',
        'TT': 'High sensitivity - significantly reduced dose needed, higher bleeding risk',
    },
    'rs8175347': {  # UGT1A1*28
        'gene': 'UGT1A1', 'drug': 'Irinotecan (chemo)',
        'evidence': 'FDA Table - severe toxicity risk',
        'TA6/TA6': 'Normal metabolism',
        'TA6/TA7': 'Intermediate - higher toxicity risk',
        'TA7/TA7': 'Poor metabolism - much higher toxicity risk',
    },
}

TIER1_CARDIOVASCULAR = {
    'rs429358': {  # APOE ε4
        'gene': 'APOE', 'condition': 'Alzheimer disease',
        'evidence': 'See "Longevity & Cognitive Aging" section for full personalized APOE interpretation',
        'TT': 'No ε4 allele - see combined interpretation with rs7412 below',
        'CT': 'ONE ε4 allele - see combined interpretation with rs7412 below',
        'CC': 'TWO ε4 copies (ε4/ε4) - very high risk - see full interpretation below',
    },
    'rs7412': {  # APOE ε2
        'gene': 'APOE', 'condition': 'Alzheimer disease',
        'evidence': 'See "Longevity & Cognitive Aging" section for full personalized APOE interpretation',
        'CC': 'No ε2 allele - see combined interpretation with rs429358 below',
        'CT': 'ONE ε2 allele - see combined interpretation with rs429358 below',
        'TT': 'TWO ε2 alleles - see combined interpretation with rs429358 below',
    },
    'rs1333049': {  # 9p21 locus
        'gene': 'CDKN2A/B', 'condition': 'Coronary artery disease',
        'evidence': 'GWAS P<10^-50 (OR=1.9, G allele OR=0.816 protective)',
        'CC': 'Elevated CAD/MI risk (~1.9x, +50% vs GG)',
        'CG': 'Moderate CAD risk (~1.5x, +25% vs GG)',
        'GG': 'Lower CAD risk (protective)',
    },
    'rs6025': {  # Factor V Leiden
        'gene': 'F5', 'condition': 'Blood clots (thrombophilia)',
        'evidence': 'ACMG Secondary Findings v3.2',
        'GG': 'Normal clotting',
        'GA': 'Factor V Leiden heterozygous (3.5-4.4x DVT risk)',
        'AA': 'Homozygous (11.4x risk) - very rare',
    },
    'rs1799963': {  # Prothrombin
        'gene': 'F2', 'condition': 'Blood clots',
        'evidence': 'ACMG Secondary Findings',
        'GG': 'Normal clotting',
        'GA': 'Carrier (2-3x clot risk)',
        'AA': 'Homozygous - very rare, very high risk',
    },
    'rs1800562': {  # HFE C282Y
        'gene': 'HFE', 'condition': 'Hemochromatosis (iron overload)',
        'evidence': 'ACMG, causes 85-90% of hemochromatosis',
        'GG': 'Normal iron regulation',
        'GA': 'Carrier (usually asymptomatic)',
        'AA': 'At risk for iron overload',
    },
    'rs1799945': {  # HFE H63D
        'gene': 'HFE', 'condition': 'Hemochromatosis (milder)',
        'evidence': 'ACMG Secondary Findings',
        'CC': 'Normal',
        'CG': 'Carrier (mild effect)',
        'GG': 'Two copies (mild iron overload risk)',
    },
}


# ============================================================================
# TIER 2: HEALTH & DISEASE RISK
# GWAS validated, replicated studies
# ============================================================================

HEALTH_METABOLIC = {
    'rs7903146': {  # TCF7L2
        'gene': 'TCF7L2', 'condition': 'Type 2 diabetes',
        'evidence': 'GWAS P<10^-100, OR 1.4-2.0',
        'CC': 'Typical risk',
        'CT': 'Elevated risk (OR ~1.4)',
        'TT': 'Higher risk (OR ~2.0)',
    },
    'rs9939609': {  # FTO
        'gene': 'FTO', 'condition': 'Obesity/BMI',
        'evidence': 'GWAS P<10^-40, +3-4kg per allele',
        'TT': 'Lower BMI tendency',
        'AT': 'Moderate BMI (+1-2kg)',
        'AA': 'Higher BMI tendency (+3-4kg)',
    },
    'rs1801282': {  # PPARG Pro12Ala
        'gene': 'PPARG', 'condition': 'Insulin sensitivity',
        'evidence': 'Complex effects on metabolism',
        'CC': 'Normal fat metabolism (Pro/Pro)',
        'CG': 'Pro/Ala - altered metabolism (conflicting evidence)',
        'GG': 'Ala/Ala - altered metabolism (conflicting evidence)',
    },
    'rs5219': {  # KCNJ11
        'gene': 'KCNJ11', 'condition': 'Type 2 diabetes',
        'evidence': 'Replicated GWAS',
        'CC': 'Typical risk',
        'CT': 'Slightly elevated',
        'TT': 'Elevated risk',
    },
    'rs10830963': {  # MTNR1B
        'gene': 'MTNR1B', 'condition': 'Fasting glucose',
        'evidence': 'Replicated GWAS',
        'CC': 'Normal',
        'CG': 'Elevated fasting glucose',
        'GG': 'Higher glucose, T2D risk',
    },
    'rs738409': {  # PNPLA3 (I148M)
        'gene': 'PNPLA3', 'condition': 'Fatty liver disease',
        'evidence': 'Strong association (NAFLD OR=3.41, Cirrhosis OR=1.86)',
        'CC': 'Normal risk (wild-type)',
        'CG': 'Elevated NAFLD risk (~2x)',
        'GG': 'High NAFLD/cirrhosis/HCC risk (I148M homozygous)',
    },
    'rs780094': {  # GCKR
        'gene': 'GCKR', 'condition': 'Triglycerides, liver fat',
        'evidence': 'GWAS validated',
        'CC': 'Normal',
        'CT': 'Moderately elevated triglycerides',
        'TT': 'Higher triglycerides, NAFLD risk',
    },
}

HEALTH_CARDIOVASCULAR_EXTENDED = {
    'rs5128': {  # APOC3
        'gene': 'APOC3', 'condition': 'Triglycerides',
        'evidence': 'Replicated studies',
        'CC': 'Typical levels',
        'CG': 'Slightly elevated',
        'GG': 'Higher triglyceride tendency',
    },
    'rs662799': {  # APOA5 (-1131T>C)
        'gene': 'APOA5', 'condition': 'Triglycerides/HDL',
        'evidence': 'Strong association (PMC8378982, MI OR=1.44)',
        'AA': 'Typical triglycerides (protective)',
        'AG': 'Elevated triglycerides (+11% per G allele)',
        'GG': 'Higher triglycerides (+36%), lower HDL, increased CVD risk',
    },
    'rs2383206': {  # CDKN2A/B
        'gene': 'CDKN2A/B', 'condition': 'Heart attack',
        'evidence': '9p21 locus',
        'AA': 'Elevated risk',
        'AG': 'Moderate risk',
        'GG': 'Lower risk',
    },
    'rs1800795': {  # IL6
        'gene': 'IL6', 'condition': 'Inflammation, CVD',
        'evidence': 'Cytokine production',
        'GG': 'Lower IL-6 production',
        'GC': 'Moderate',
        'CC': 'Higher IL-6, inflammation',
    },
    # Note: rs4420638 (APOE region) moved to LONGEVITY_COGNITIVE_AGING section for consolidated interpretation
}

HEALTH_CANCER = {
    'rs1042522': {  # TP53
        'gene': 'TP53', 'condition': 'Tumor suppressor (Arg72Pro)',
        'evidence': 'Functional variant',
        'GG': 'Arg/Arg - common variant, slightly shorter lifespan',
        'CG': 'Arg/Pro - intermediate',
        'CC': 'Pro/Pro - may live ~3 years longer, better chemo response',
    },
    'rs17849079': {  # PIK3CA
        'gene': 'PIK3CA', 'condition': 'Cowden syndrome',
        'evidence': 'PI3K-AKT pathway',
        'GG': 'Normal risk',
        'GT': 'Slightly elevated risk',
        'TT': 'Elevated risk',
    },
    'rs2227983': {  # EGFR
        'gene': 'EGFR', 'condition': 'Lung cancer (protective)',
        'evidence': 'Growth factor receptor',
        'AA': 'Higher lung cancer risk if smoker',
        'AG': 'Moderate protection',
        'GG': 'Better protection against lung cancer',
    },
    'rs4073': {  # CXCL8
        'gene': 'CXCL8', 'condition': 'Various cancers',
        'evidence': 'Inflammation marker',
        'AA': 'Lower IL-8',
        'AT': 'Moderate',
        'TT': 'Higher IL-8, increased cancer risk',
    },
    'rs1800896': {  # IL10
        'gene': 'IL10', 'condition': 'Immune function',
        'evidence': 'Anti-inflammatory cytokine',
        'AA': 'Higher IL-10 (protective)',
        'AG': 'Moderate',
        'GG': 'Lower IL-10',
    },
    'rs4986790': {  # TLR4
        'gene': 'TLR4', 'condition': 'Immune response, infection risk',
        'evidence': 'Innate immunity',
        'AA': 'Normal response',
        'AG': 'Altered immune response',
        'GG': 'Different response pattern',
    },
}

HEALTH_COVID19 = {
    'rs10490770': {  # COVID-19 severity (3p21.31)
        'gene': 'LZTFL1', 'condition': 'COVID-19 severity',
        'evidence': 'GWAS - Nature 2020',
        'AA': 'Lower risk of severe COVID',
        'AG': 'Moderate risk (~1.7x)',
        'GG': '2x risk of severe COVID-19',
    },
    'rs657152': {  # ABO blood type
        'gene': 'ABO', 'condition': 'Blood type (COVID-19 susceptibility)',
        'evidence': 'Blood type O protective',
        'AA': 'Blood type A (higher COVID risk)',
        'AC': 'Blood type A or AB',
        'CC': 'Blood type O (lower COVID risk)',
    },
}

HEALTH_AUTOIMMUNE = {
    'rs2476601': {  # PTPN22
        'gene': 'PTPN22', 'condition': 'Autoimmune (T1D, RA, etc)',
        'evidence': 'Strong autoimmune association',
        'GG': 'Typical risk',
        'AG': 'Elevated autoimmune risk',
        'AA': 'Higher autoimmune risk',
    },
    'rs1800629': {  # TNF
        'gene': 'TNF', 'condition': 'TNF-alpha production',
        'evidence': 'Pro-inflammatory cytokine',
        'GG': 'Normal TNF levels',
        'GA': 'Higher TNF production',
        'AA': 'Much higher TNF (pro-inflammatory)',
    },
    'rs2187668': {  # HLA-DQ
        'gene': 'HLA-DQA1', 'condition': 'Celiac disease',
        'evidence': 'HLA-DQ region',
        'CC': 'Very low risk (<1%)',
        'CT': 'Low-moderate risk',
        'TT': 'Elevated risk',
    },
    'rs6457620': {  # HLA region
        'gene': 'HLA region', 'condition': 'Rheumatoid arthritis',
        'evidence': 'HLA association',
        'CC': 'Lower risk',
        'CT': 'Moderate risk',
        'TT': 'Higher RA risk',
    },
}

HEALTH_BONE_KIDNEY = {
    'rs2234693': {  # ESR1
        'gene': 'ESR1', 'condition': 'Osteoporosis',
        'evidence': 'Estrogen receptor',
        'CC': 'Lower bone density tendency',
        'CT': 'Moderate',
        'TT': 'Better bone density',
    },
    'rs4293393': {  # UMOD
        'gene': 'UMOD', 'condition': 'Chronic kidney disease',
        'evidence': 'GWAS validated (PNAS 2022, CKD OR=1.25, T allele=risk)',
        'AA': 'Elevated CKD risk, salt-sensitive hypertension (T/T on research strand)',
        'AG': 'Moderate CKD risk',
        'GG': 'Lower CKD risk (protective, C/C on research strand, lower uromodulin)',
    },
    'rs1799983': {  # NOS3
        'gene': 'NOS3', 'condition': 'Nitric oxide, blood pressure',
        'evidence': 'Endothelial function',
        'GG': 'Normal NO production',
        'GT': 'Slightly reduced',
        'TT': 'Lower NO, hypertension risk',
    },
}

# ============================================================================
# LONGEVITY & COGNITIVE AGING - FOXO3 + APOE
# Most replicated longevity genes - interpret together
# ============================================================================

LONGEVITY_COGNITIVE_AGING = {
    # FOXO3 - The Longevity Gene (most replicated across all populations)
    'rs2802292': {  # FOXO3
        'gene': 'FOXO3', 'condition': 'Longevity, healthy aging',
        'evidence': 'Most replicated longevity gene worldwide - centenarian association across all populations',
        'TT': '🌟 STRONG longevity (~1.8x centenarian odds) - stress resistance, healthy aging',
        'GT': '✓ MODERATE longevity (~1.3x centenarian odds) - good aging trajectory',
        'GG': 'Typical lifespan',
    },
    'rs2764264': {  # FOXO3
        'gene': 'FOXO3', 'condition': 'Longevity (male-specific effect)',
        'evidence': 'Centenarian studies - stronger effect in males',
        'TT': '🌟 STRONG longevity association',
        'CT': '✓ MODERATE longevity association',
        'CC': 'Typical lifespan',
    },

    # APOE - Alzheimer's & Cognitive Aging (interpret rs429358 + rs7412 together!)
    'rs429358': {  # APOE ε4 - THE CRITICAL ONE
        'gene': 'APOE', 'condition': 'Alzheimer disease, cognitive aging, lifespan',
        'evidence': 'Must combine with rs7412 to determine APOE type - see personalized interpretation below',
        'TT': 'No ε4 allele (need rs7412: if TT=ε2/ε2 [0.6x risk], if CT=ε2/ε3 [0.6x risk], if CC=ε3/ε3 [1.0x risk])',
        'CT': 'ONE ε4 allele (need rs7412: if CT=ε2/ε4 [~2.6x risk], if CC=ε3/ε4 [~3.2x risk])',
        'CC': '⚠⚠ TWO ε4 copies = ε4/ε4 - very high Alzheimer risk (~12x, 61x early-onset), very rare (~2%)',
    },
    'rs7412': {  # APOE ε2 - THE PROTECTIVE ONE
        'gene': 'APOE', 'condition': 'Alzheimer protection, longevity',
        'evidence': 'Must combine with rs429358 to determine APOE type - see personalized interpretation below',
        'TT': 'TWO ε2 alleles (need rs429358: if TT=ε2/ε2 [0.6x risk, longevity])',
        'CT': 'ONE ε2 allele (need rs429358: if TT=ε2/ε3 [0.6x risk], if CT=ε2/ε4 [~2.6x risk])',
        'CC': 'No ε2 allele (need rs429358: if TT=ε3/ε3 [1.0x risk], if CT=ε3/ε4 [~3.2x risk], if CC=ε4/ε4 [~12x risk])',
    },

    # APOE region proxy marker
    'rs4420638': {  # APOE region
        'gene': 'APOE region', 'condition': 'Alzheimer risk (additional marker), cholesterol',
        'evidence': 'APOE region proxy - correlates with ε4 status',
        'GG': '⚠ Higher Alzheimer risk (2x+), higher LDL cholesterol',
        'AG': '⚠ Elevated Alzheimer risk (~3x), higher heart disease risk (1.4x)',
        'AA': '✓ Normal/average Alzheimer risk',
    },

    # Other longevity markers
    'rs3764261': {  # CETP
        'gene': 'CETP', 'condition': 'HDL cholesterol, longevity',
        'evidence': 'Higher HDL, longevity association (PMC3293889)',
        'AA': '✓ Lower CETP activity - higher HDL ("good" cholesterol), longevity association',
        'CA': 'Reduced CETP activity - higher HDL',
        'CC': 'Typical CETP activity - typical HDL levels',
    },
    'rs9536314': {  # KLOTHO (F352V, KL-VS variant)
        'gene': 'KLOTHO', 'condition': 'Anti-aging, cognition',
        'evidence': 'Anti-aging gene - klotho protein (PMC4978356)',
        'TT': 'Typical klotho function',
        'GT': '✓ KL-VS heterozygote - better cognition, may enhance longevity',
        'GG': 'KL-VS homozygote - may reduce lifespan (rare)',
    },
}

# ============================================================================
# LD-BASED SNP PREDICTION
# ============================================================================

def predict_foxo3_from_ld(snps_obj):
    """
    Predict rs2802292 (FOXO3 longevity SNP) from proxy SNPs in high LD

    Based on research: rs2802292 has 7+ proxy SNPs with r² 0.89-0.99
    These SNPs form a longevity haplotype that's highly conserved.
    """
    # Proxy SNPs in strong LD with rs2802292 (r² > 0.89)
    # Source: PMC6606898, various FOXO3 longevity studies
    FOXO3_PROXIES = {
        'rs2802288': {  # Very strong proxy (r² ≈ 0.95-0.99)
            'TT': 'TT',  # If proxy is TT, target is TT (longevity allele)
            'GT': 'GT',  # Heterozygous
            'GG': 'GG',  # Non-longevity allele
            'r2': 0.97
        },
        'rs2764264': {  # Another FOXO3 longevity SNP (moderate LD, r² ≈ 0.7-0.8)
            'TT': 'TT',  # Both are longevity markers
            'CT': 'GT',  # Heterozygous
            'CC': 'GG',
            'r2': 0.75
        },
        'rs12202234': {  # Strong LD (r² ≈ 0.90-0.95)
            'TT': 'TT',
            'CT': 'GT',
            'CC': 'GG',
            'r2': 0.92
        },
        'rs3800230': {  # Strong LD (r² ≈ 0.89-0.93)
            'AA': 'TT',  # Different alleles, but correlated
            'AG': 'GT',
            'GG': 'GG',
            'r2': 0.91
        },
    }

    predictions = []

    for proxy_rsid, proxy_info in FOXO3_PROXIES.items():
        proxy_data = snps_obj.snps[snps_obj.snps.index == proxy_rsid]

        if len(proxy_data) > 0 and not proxy_data['genotype'].isna().all():
            proxy_gt = str(proxy_data['genotype'].iloc[0])

            # Normalize (handle reversed genotypes)
            if proxy_gt in ['TC', 'CT']:
                proxy_gt = 'CT' if 'C' in proxy_gt and 'T' in proxy_gt else proxy_gt
            elif proxy_gt in ['GA', 'AG']:
                proxy_gt = 'AG' if 'A' in proxy_gt and 'G' in proxy_gt else proxy_gt

            if proxy_gt in proxy_info:
                predicted_gt = proxy_info[proxy_gt]
                r2 = proxy_info['r2']
                predictions.append({
                    'genotype': predicted_gt,
                    'r2': r2,
                    'proxy_rsid': proxy_rsid,
                    'proxy_genotype': proxy_gt
                })

    if not predictions:
        return None

    # Use the prediction from the highest r² proxy
    best_prediction = max(predictions, key=lambda x: x['r2'])

    # Check for consensus (all predictions agree)
    consensus = all(p['genotype'] == best_prediction['genotype'] for p in predictions)

    return {
        'genotype': best_prediction['genotype'],
        'r2': best_prediction['r2'],
        'num_proxies': len(predictions),
        'consensus': consensus,
        'method': 'LD-based prediction',
        'proxies': predictions
    }

# ============================================================================
# PERSONALIZED FOXO3 + APOE INTERPRETATION
# ============================================================================

def interpret_apoe_foxo3(snps_obj):
    """
    Interpret user's specific APOE and FOXO3 genotypes with personalized recommendations
    """
    # Get APOE genotypes
    rs429358_data = snps_obj.snps[snps_obj.snps.index == 'rs429358']
    rs7412_data = snps_obj.snps[snps_obj.snps.index == 'rs7412']

    rs429358_gt = None
    rs7412_gt = None

    if len(rs429358_data) > 0 and not rs429358_data['genotype'].isna().all():
        rs429358_gt = str(rs429358_data['genotype'].iloc[0])
    if len(rs7412_data) > 0 and not rs7412_data['genotype'].isna().all():
        rs7412_gt = str(rs7412_data['genotype'].iloc[0])

    # Get FOXO3 genotypes (with LD-based prediction fallback)
    rs2802292_data = snps_obj.snps[snps_obj.snps.index == 'rs2802292']
    rs2764264_data = snps_obj.snps[snps_obj.snps.index == 'rs2764264']

    rs2802292_gt = None
    rs2802292_predicted = False
    rs2802292_prediction_info = None

    # Try direct genotyping first
    if len(rs2802292_data) > 0 and not rs2802292_data['genotype'].isna().all():
        rs2802292_gt = str(rs2802292_data['genotype'].iloc[0])
    else:
        # Try LD-based prediction
        prediction = predict_foxo3_from_ld(snps_obj)
        if prediction:
            rs2802292_gt = prediction['genotype']
            rs2802292_predicted = True
            rs2802292_prediction_info = prediction

    rs2764264_gt = None
    if len(rs2764264_data) > 0 and not rs2764264_data['genotype'].isna().all():
        rs2764264_gt = str(rs2764264_data['genotype'].iloc[0])

    # Determine APOE type
    apoe_type = None
    apoe_risk = None
    apoe_description = None
    apoe_recommendations = []

    if rs429358_gt and rs7412_gt:
        # Normalize genotypes (handle reversed)
        if rs429358_gt == 'TC': rs429358_gt = 'CT'
        if rs7412_gt == 'TC': rs7412_gt = 'CT'

        # Map to APOE type
        apoe_map = {
            ('TT', 'TT'): ('ε2/ε2', '0.6x', '🌟 PROTECTIVE', 'Very rare (~1%). Strong Alzheimer protection and enriched in centenarians. However, can cause type III hyperlipoproteinemia - monitor triglycerides.'),
            ('TT', 'CT'): ('ε2/ε3', '0.6x', '✓ PROTECTIVE', 'Protective variant (~10-15% of population). Lower Alzheimer risk than average.'),
            ('CT', 'CT'): ('ε2/ε4', '~2.6x', '⚠ INTERMEDIATE', 'Rare combination (~2%). You have both the protective ε2 and risk ε4 alleles. The ε2 partially offsets ε4 risk, resulting in intermediate Alzheimer risk.'),
            ('TT', 'CC'): ('ε3/ε3', '1.0x', '✓ BASELINE', 'Most common type (~60-70% of population). Normal/average Alzheimer risk.'),
            ('CT', 'CC'): ('ε3/ε4', '~3.2x', '⚠ ELEVATED', 'Common variant (~25% of population). 3x higher Alzheimer risk, 1.4x heart disease risk. Focus on prevention!'),
            ('CC', 'CC'): ('ε4/ε4', '~12x', '⚠⚠ HIGH', 'Very rare (~2%). Highest Alzheimer risk (12x late-onset, 61x early-onset). Strong prevention strategy recommended.'),
        }

        key = (rs429358_gt, rs7412_gt)
        if key in apoe_map:
            apoe_type, apoe_risk, apoe_icon, apoe_description = apoe_map[key]

            # Recommendations based on APOE type
            if 'ε4' in apoe_type:
                apoe_recommendations = [
                    "🏃 Exercise regularly (150+ min/week aerobic + strength training)",
                    "🥗 Mediterranean diet (high fish, olive oil, vegetables, low red meat)",
                    "💤 Prioritize sleep (7-9 hours, treat sleep apnea if present)",
                    "🧠 Cognitive training and lifelong learning",
                    "❤️ Control cardiovascular risk factors (BP, cholesterol, diabetes)",
                    "🚭 Avoid smoking and excessive alcohol",
                    "👥 Stay socially engaged",
                ]
                if apoe_type == 'ε4/ε4':
                    apoe_recommendations.insert(0, "⚠️ Consider genetic counseling and early/frequent cognitive screening")
            elif 'ε2' in apoe_type:
                apoe_recommendations = [
                    "🌟 You have genetic protection against Alzheimer's!",
                    "📊 Monitor triglycerides (ε2 can rarely cause elevated levels)",
                    "💪 Continue healthy lifestyle to maximize longevity benefits",
                ]
            else:  # ε3/ε3
                apoe_recommendations = [
                    "✓ Average Alzheimer risk - standard prevention applies",
                    "🏃 Regular exercise and healthy diet still beneficial",
                    "🧠 Cognitive engagement supports healthy aging",
                ]

    # Interpret FOXO3
    foxo3_longevity = None
    foxo3_description = None

    if rs2802292_gt:
        if rs2802292_gt in ['TT']:
            foxo3_longevity = '🌟 STRONG'
            foxo3_description = 'EXCELLENT longevity genetics! ~1.8x odds of reaching 100. Associated with stress resistance and healthy aging.'
        elif rs2802292_gt in ['GT', 'TG']:
            foxo3_longevity = '✓ MODERATE'
            foxo3_description = 'Good longevity genetics. ~1.3x odds of reaching 100. Favorable aging trajectory.'
        else:  # GG
            foxo3_longevity = 'TYPICAL'
            foxo3_description = 'Average longevity - lifestyle factors still matter greatly!'

    # Only print if we have SOMETHING to show
    if not (foxo3_longevity or apoe_type or rs429358_gt or rs7412_gt):
        return  # No data available, skip this section

    # Print interpretation
    print(f"\n{Colors.BOLD}{Colors.CYAN}{'='*80}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}  YOUR PERSONALIZED LONGEVITY & COGNITIVE AGING PROFILE{Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}{'='*80}{Colors.END}\n")

    # FOXO3 Section
    if foxo3_longevity:
        if rs2802292_predicted:
            # Show predicted result with details
            print(f"{Colors.BOLD}🧬 FOXO3 Longevity Gene (rs2802292: {rs2802292_gt}*):{Colors.END}")
            print(f"  {foxo3_description}\n")

            # Show prediction details
            info = rs2802292_prediction_info
            confidence_color = Colors.GREEN if info['r2'] > 0.9 else Colors.YELLOW if info['r2'] > 0.8 else Colors.RED

            print(f"  {Colors.CYAN}* Predicted from {info['num_proxies']} proxy SNP(s) (not directly genotyped){Colors.END}")
            print(f"  {Colors.GRAY}Method: {info['method']}{Colors.END}")
            print(f"  {Colors.GRAY}Confidence: {confidence_color}r² = {info['r2']:.2f}{Colors.END} {Colors.GRAY}(LD correlation){Colors.END}")

            if info['consensus']:
                print(f"  {Colors.GREEN}✓ All proxy SNPs agree on this prediction{Colors.END}")
            else:
                print(f"  {Colors.YELLOW}⚠ Multiple predictions - using highest r² proxy{Colors.END}")

            # Show which proxy was used
            best = max(info['proxies'], key=lambda x: x['r2'])
            print(f"  {Colors.GRAY}Primary proxy: {best['proxy_rsid']} ({best['proxy_genotype']}) → r² = {best['r2']:.2f}{Colors.END}")

            if info['r2'] < 0.85:
                print(f"\n  {Colors.YELLOW}⚠ Note: Moderate confidence prediction - for curiosity only{Colors.END}")
                print(f"  {Colors.YELLOW}Consider 23andMe or whole genome sequencing for definitive results{Colors.END}")

            print()
        else:
            # Directly genotyped
            print(f"{Colors.BOLD}🧬 FOXO3 Longevity Gene (rs2802292: {rs2802292_gt}):{Colors.END}")
            print(f"  {foxo3_description}")
            print(f"  {Colors.GREEN}✓ Directly genotyped (high confidence){Colors.END}\n")
    else:
        print(f"{Colors.BOLD}🧬 FOXO3 Longevity Gene:{Colors.END}")
        print(f"  {Colors.GRAY}⚠️ Data not available - rs2802292 not found and no proxy SNPs detected{Colors.END}")
        print(f"  {Colors.GRAY}(Most replicated longevity gene - associated with centenarian status){Colors.END}")
        print(f"\n  {Colors.GRAY}Missing proxies checked: rs2802288, rs2764264, rs12202234, rs3800230{Colors.END}")
        print(f"  {Colors.GRAY}Recommendation: Consider 23andMe or whole genome sequencing{Colors.END}\n")

    # APOE Section - Complete data
    if apoe_type:
        risk_color = Colors.GREEN if '0.6x' in apoe_risk or '1.0x' in apoe_risk else Colors.YELLOW if '2.6x' in apoe_risk else Colors.RED

        print(f"{Colors.BOLD}🧠 APOE Type (rs429358={rs429358_gt} + rs7412={rs7412_gt}):{Colors.END}")
        print(f"  {Colors.BOLD}Type:{Colors.END} {Colors.BOLD}{apoe_type}{Colors.END}")
        print(f"  {Colors.BOLD}Alzheimer Risk:{Colors.END} {risk_color}{apoe_risk}{Colors.END} vs. average")
        print(f"  {Colors.BOLD}Assessment:{Colors.END} {apoe_description}\n")

        if apoe_recommendations:
            print(f"{Colors.BOLD}📋 Personalized Recommendations:{Colors.END}")
            for rec in apoe_recommendations:
                print(f"  {rec}")
            print()

    # APOE Section - Incomplete data
    elif rs429358_gt or rs7412_gt:
        print(f"{Colors.BOLD}🧠 APOE Genotype (Incomplete Data):{Colors.END}")
        if rs429358_gt and not rs7412_gt:
            print(f"  {Colors.YELLOW}⚠ Only rs429358 available ({rs429358_gt}) - missing rs7412{Colors.END}")
            print(f"  Cannot determine complete APOE type without both SNPs.")
            if rs429358_gt == 'CC':
                print(f"  {Colors.RED}Note: CC indicates ε4/ε4 (highest risk) - consider follow-up testing{Colors.END}")
            elif rs429358_gt in ['CT', 'TC']:
                print(f"  Note: You have one ε4 allele - could be ε2/ε4 or ε3/ε4 depending on rs7412")
            else:  # TT
                print(f"  Note: No ε4 allele - could be ε2/ε2, ε2/ε3, or ε3/ε3 depending on rs7412")
        elif rs7412_gt and not rs429358_gt:
            print(f"  {Colors.YELLOW}⚠ Only rs7412 available ({rs7412_gt}) - missing rs429358{Colors.END}")
            print(f"  Cannot determine complete APOE type without both SNPs.")
            if rs7412_gt == 'TT':
                print(f"  Note: You have TWO ε2 alleles - likely protective if rs429358 confirms")
            elif rs7412_gt in ['CT', 'TC']:
                print(f"  Note: You have ONE ε2 allele - could be ε2/ε3 or ε2/ε4 depending on rs429358")
            else:  # CC
                print(f"  Note: No ε2 allele - could be ε3/ε3, ε3/ε4, or ε4/ε4 depending on rs429358")
        print(f"\n  {Colors.GRAY}Recommendation: Consider 23andMe, AncestryDNA, or clinical APOE testing{Colors.END}")
        print(f"  {Colors.GRAY}for complete genotype information.{Colors.END}\n")

    # APOE Section - No data at all
    else:
        print(f"{Colors.BOLD}🧠 APOE Genotype:{Colors.END}")
        print(f"  {Colors.GRAY}⚠️ Data not available - neither rs429358 nor rs7412 found in your DNA file{Colors.END}")
        print(f"  {Colors.GRAY}APOE is the strongest genetic risk factor for Alzheimer's disease.{Colors.END}")
        print(f"\n  {Colors.GRAY}Recommendation: Consider 23andMe, AncestryDNA, or clinical APOE testing.{Colors.END}")
        print(f"  {Colors.GRAY}APOE testing is particularly important if you have family history of Alzheimer's.{Colors.END}\n")

    # Combined interpretation
    if foxo3_longevity and apoe_type:
        print(f"{Colors.BOLD}🔬 Combined Profile:{Colors.END}")

        if 'STRONG' in foxo3_longevity and 'ε4' in apoe_type:
            print(f"  {Colors.YELLOW}You have FOXO3 longevity genetics which may partially offset APOE ε4 risk.{Colors.END}")
            print(f"  {Colors.YELLOW}Focus on prevention strategies - your longevity genes give you time to benefit!{Colors.END}")
        elif 'STRONG' in foxo3_longevity:
            print(f"  {Colors.GREEN}Excellent combination! Strong longevity genetics + favorable/normal cognitive aging profile.{Colors.END}")
        elif 'ε2' in apoe_type:
            print(f"  {Colors.GREEN}Great cognitive aging protection from APOE ε2!{Colors.END}")

        print()

    print(f"{Colors.GRAY}{'─'*80}{Colors.END}\n")

HEALTH_LONGEVITY = LONGEVITY_COGNITIVE_AGING  # Keep alias for backward compatibility

HEALTH_COGNITIVE = {
    # Note: APOE is now in LONGEVITY_COGNITIVE_AGING section above for consolidated interpretation

    # Working Memory & Executive Function
    'rs4680_cognition': {  # COMT Val158Met
        'gene': 'COMT', 'condition': 'Working memory, executive function, attention',
        'evidence': 'Val158Met - most studied cognitive SNP, validated across populations',
        'AA': 'Met/Met - "Worrier" - better memory, attention, information processing under normal conditions',
        'GA': 'Val/Met - intermediate cognitive performance',
        'GG': 'Val/Val - "Warrior" - modest reduction in executive function, better under stress',
    },
    'rs1006737_cognition': {  # CACNA1C
        'gene': 'CACNA1C', 'condition': 'Working memory, prefrontal function',
        'evidence': 'L-type calcium channel - validated in multiple studies',
        'GG': 'Better working memory performance',
        'AG': 'Intermediate working memory',
        'AA': 'Impaired working memory, reduced prefrontal efficiency',
    },
    'rs1044396': {  # CHRNA4
        'gene': 'CHRNA4', 'condition': 'Attention, working memory, parietal cortex function',
        'evidence': 'Nicotinic receptor - validated for attention and cognitive control',
        'CC': 'Reduced parietal cortex activity during attention tasks',
        'CT': 'Intermediate attention performance',
        'TT': 'Robust parietal cortex activity, enhanced attention and working memory',
    },

    # Memory Performance
    'rs17070145': {  # KIBRA
        'gene': 'KIBRA', 'condition': 'Episodic memory performance',
        'evidence': 'Replicated association with memory recall',
        'CC': 'Reduced memory abilities',
        'CT': 'Increased memory performance (T allele +24% recall)',
        'TT': 'Greatly increased memory performance',
    },
    'rs6265_cognition': {  # BDNF
        'gene': 'BDNF', 'condition': 'Learning, memory consolidation',
        'evidence': 'Val66Met affects learning and neuroplasticity',
        'GG': 'Val/Val - better learning, memory, neuroplasticity',
        'GA': 'Val/Met - impaired motor learning, reduced BDNF secretion',
        'AA': 'Met/Met - reduced learning efficiency, faster cognitive decline with age',
    },

    # Cognitive Decline & Alzheimer's Risk
    'rs11136000': {  # CLU (Clusterin)
        'gene': 'CLU', 'condition': 'Cognitive decline, memory, Alzheimer\'s disease',
        'evidence': '2nd strongest AD risk SNP after APOE - validated in multiple GWAS',
        'TT': 'Lower risk of cognitive decline, better memory maintenance',
        'CT': 'Intermediate risk - faster memory decline in presymptomatic stages',
        'CC': 'Higher risk - accelerated cognitive decline, poorer memory scores, faster MCI→AD progression',
    },
    'rs9331896': {  # CLU (Clusterin)
        'gene': 'CLU', 'condition': 'Cognitive decline, hippocampal function',
        'evidence': 'Strongest CLU SNP for AD risk - affects CLU expression in hippocampus',
        'CC': 'Lower Alzheimer\'s risk, better memory maintenance (protective)',
        'CT': 'Intermediate risk of cognitive decline (OR 1.18)',
        'TT': 'Higher AD risk, faster presymptomatic memory decline, altered hippocampal CLU expression',
    },

    # Other markers
    'rs363050': {  # SNAP25
        'gene': 'SNAP25', 'condition': 'Intelligence, IQ',
        'evidence': 'Original findings failed to replicate - interpret with caution',
        'AA': 'Original (unreplicated) study suggested +2.8 PIQ vs GG',
        'AG': 'Intermediate (unreplicated findings)',
        'GG': 'Reference genotype (note: original IQ findings did not replicate)',
    },
}

HEALTH_PAIN = {
    'rs6746030': {  # SCN9A
        'gene': 'SCN9A', 'condition': 'Pain sensitivity',
        'evidence': 'Sodium channel - pain perception',
        'GG': 'Normal pain sensitivity',
        'GA': 'Reduced pain sensitivity',
        'AA': 'Lower pain sensitivity',
    },
    'rs6267': {  # COMT
        'gene': 'COMT', 'condition': 'Schizophrenia associations (research mixed/inconclusive)',
        'evidence': 'COMT variant - some studies suggest schizophrenia association, others find no link',
        'GG': 'Common genotype (~83% frequency, research on schizophrenia association is inconclusive)',
        'GT': 'Less common (~16%, TT genotype codes for similar enzyme activity as GT)',
        'TT': 'Rare genotype (~1%, research on schizophrenia association is inconclusive)',
    },
    'rs1799971_pain': {  # OPRM1
        'gene': 'OPRM1', 'condition': 'Pain, opioid response',
        'evidence': 'Opioid receptor efficacy',
        'AA': 'Normal pain response, better opioid effect',
        'AG': 'Intermediate',
        'GG': 'Higher pain sensitivity, reduced opioid efficacy',
    },
}

HEALTH_ADDICTION = {
    'rs1800497_addiction': {  # DRD2
        'gene': 'DRD2', 'condition': 'Addiction susceptibility',
        'evidence': 'Reward deficiency syndrome',
        'CC': 'Normal D2 receptors',
        'CT': 'Reduced receptors - addiction risk',
        'TT': 'Low D2 - higher addiction/substance abuse risk',
    },
    'rs279858': {  # GABRA2
        'gene': 'GABRA2', 'condition': 'Alcohol dependence',
        'evidence': 'Alcoholism association',
        'AA': 'Lower alcoholism risk',
        'AG': 'Elevated risk (slower alcohol response)',
        'GG': 'Higher alcoholism risk (slower alcohol response)',
    },
    'rs2023239': {  # CNR1
        'gene': 'CNR1', 'condition': 'Cannabis dependence',
        'evidence': 'Cannabinoid receptor',
        'CC': 'Typical',
        'CT': 'Moderate risk',
        'TT': 'Higher cannabis dependence risk',
    },
    'rs806378': {  # CNR1
        'gene': 'CNR1', 'condition': 'Substance dependence',
        'evidence': 'Addiction studies',
        'CC': 'Lower risk',
        'CT': 'Moderate',
        'TT': 'Elevated addiction risk',
    },
    'rs324420_addiction': {  # FAAH
        'gene': 'FAAH', 'condition': 'Substance use disorder',
        'evidence': 'Endocannabinoid system',
        'CC': 'Normal risk',
        'CA': 'Intermediate',
        'AA': 'Significantly increased substance use disorder risk',
    },
}


# ============================================================================
# MOOD, MENTAL HEALTH & NEUROTRANSMITTERS
# Comprehensive collection
# ============================================================================

MOOD_SEROTONIN = {
    'rs6295': {  # HTR1A
        'gene': 'HTR1A', 'system': 'Serotonin 1A receptor',
        'effect': 'Depression, anxiety, SSRI response',
        'CC': 'Normal receptor function',
        'CG': 'Altered SSRI response',
        'GG': 'Lower expression, depression/anxiety risk',
    },
    'rs6313': {  # HTR2A
        'gene': 'HTR2A', 'system': 'Serotonin 2A receptor',
        'effect': 'SSRI response, OCD, psychedelics',
        'GG': 'Better SSRI response',
        'GA': 'Intermediate',
        'AA': 'Poorer SSRI response',
    },
    'rs7997012': {  # HTR2A
        'gene': 'HTR2A', 'system': 'Serotonin 2A',
        'effect': 'Depression treatment response (citalopram)',
        'AA': '18% better response to citalopram',
        'AG': 'Normal response',
        'GG': '18% worse response to citalopram',
    },
    'rs6311': {  # HTR2A promoter
        'gene': 'HTR2A', 'system': 'Serotonin 2A expression',
        'effect': 'SSRI sexual dysfunction, suicide risk',
        'TT': 'Normal (lower) risk of SSRI sexual dysfunction, protective against suicide',
        'CT': 'Normal risk',
        'CC': '3.6x risk of sexual dysfunction on SSRIs, higher suicide risk',
    },
    'rs4570625': {  # TPH2
        'gene': 'TPH2', 'system': 'Serotonin synthesis',
        'effect': 'Anxiety-related traits, placebo response',
        'TT': 'Normal',
        'GT': 'Normal',
        'GG': 'Higher anxiety-related traits, greater placebo response',
    },
    'rs25531': {  # SLC6A4 (5-HTTLPR)
        'gene': 'SLC6A4', 'system': 'Serotonin transporter',
        'effect': 'Stress resilience, life satisfaction',
        'GG': 'More resilient to stress, optimistic, higher life satisfaction',
        'AG': 'Intermediate',
        'AA': 'Lower serotonin, slightly less happy, may need more support',
    },
    'rs25532': {  # SLC6A4
        'gene': 'SLC6A4', 'system': 'Serotonin transporter',
        'effect': 'Regulatory region, OCD association',
        'TT': 'Normal',
        'CT': 'May be part of OCD haplotype',
        'CC': 'Higher expressing allele, may be part of OCD haplotype',
    },
}

MOOD_DOPAMINE = {
    'rs4680': {  # COMT Val158Met
        'gene': 'COMT', 'system': 'Dopamine clearance',
        'effect': 'Warrior vs Worrier, stress response',
        'GG': 'Val/Val - fast clearance (Warrior, stress-resilient)',
        'GA': 'Val/Met - intermediate',
        'AA': 'Met/Met - slow clearance (Worrier, anxious, better cognition)',
    },
    'rs6269': {  # COMT
        'gene': 'COMT', 'system': 'Dopamine metabolism',
        'effect': 'Pain sensitivity, anxiety',
        'AA': 'Lower COMT activity',
        'AG': 'Intermediate',
        'GG': 'Higher activity',
    },
    'rs4633': {  # COMT
        'gene': 'COMT', 'system': 'Dopamine',
        'effect': 'Endometrial cancer risk (affects COMT expression)',
        'CC': 'Normal cancer risk (magnitude 0)',
        'CT': 'Higher endometrial cancer risk (magnitude 2)',
        'TT': 'Higher endometrial cancer risk (magnitude 2, OR 2.39)',
    },
    'rs165599': {  # COMT region
        'gene': 'COMT', 'system': 'Dopamine',
        'effect': 'Bipolar, schizophrenia',
        'GG': 'Elevated risk (magnitude 1.5)',
        'GA': 'Moderate risk (magnitude 1)',
        'AA': 'Lower risk (magnitude 0)',
    },
    'rs1800497': {  # ANKK1/DRD2 (Taq1A)
        'gene': 'DRD2', 'system': 'Dopamine D2 receptor',
        'effect': 'Reward sensitivity, addiction',
        'CC': 'Normal D2 density',
        'CT': 'Reduced receptor (reward deficiency)',
        'TT': 'Lower D2 (addiction/impulse risk)',
    },
    'rs6277': {  # DRD2
        'gene': 'DRD2', 'system': 'Dopamine D2 receptor',
        'effect': 'Schizophrenia, addiction',
        'CC': 'Higher D2 expression, 1.6x schizophrenia risk',
        'CT': 'Intermediate D2 expression, 1.4x schizophrenia risk',
        'TT': 'Lower D2 expression (T allele decreases mRNA), normal schizophrenia risk',
    },
    'rs1076560': {  # DRD2
        'gene': 'DRD2', 'system': 'Dopamine D2',
        'effect': 'Working memory, alcoholism risk',
        'CC': 'Typical',
        'AC': 'Moderate (1.3x alcoholism risk)',
        'AA': 'Altered working memory',
    },
    'rs1800955': {  # DRD4
        'gene': 'DRD4', 'system': 'Dopamine D4 receptor',
        'effect': 'Novelty-seeking, ADHD',
        'TT': 'Normal',
        'CT': 'Increased novelty-seeking susceptibility',
        'CC': 'Increased novelty-seeking susceptibility',
    },
    'rs936461': {  # DRD4
        'gene': 'DRD4', 'system': 'Dopamine D4',
        'effect': 'ADHD, personality',
        'GG': 'Typical',
        'GA': 'Moderate',
        'AA': 'ADHD association',
    },
    'rs1611115': {  # DBH
        'gene': 'DBH', 'system': 'Dopamine to norepinephrine',
        'effect': 'ADHD, mood, blood pressure',
        'TT': 'Lowest DBH enzyme activity (~4.1), ADHD/impulsiveness risk',
        'CT': 'Intermediate activity (~25.2)',
        'CC': 'Highest DBH enzyme activity (~48.1)',
    },
    'rs1108580': {  # DBH
        'gene': 'DBH', 'system': 'Dopamine beta-hydroxylase',
        'effect': 'ADHD, Parkinson, cocaine dependence associations',
        'GG': 'Variant genotype (non-pathogenic)',
        'GA': 'Heterozygous',
        'AA': 'Reference/common genotype',
    },
}

MOOD_NEUROPLASTICITY = {
    'rs6265': {  # BDNF Val66Met
        'gene': 'BDNF', 'system': 'Brain-derived neurotrophic factor',
        'effect': 'Learning, memory, stress response',
        'GG': 'Val/Val - better memory, learning',
        'GA': 'Val/Met - intermediate, impaired motor learning',
        'AA': 'Met/Met - reduced BDNF, impaired memory but may be depression resistant',
    },
    'rs11030104': {  # BDNF
        'gene': 'BDNF', 'system': 'Neuroplasticity',
        'effect': 'BMI, Alzheimer disease, bipolar disorder',
        'AA': 'Higher BMI tendency (A is risk allele)',
        'AG': 'Moderate BMI increase, may affect Alzheimer (non-ApoE4)',
        'GG': 'Lower BMI tendency',
    },
    'rs908867': {  # BDNF
        'gene': 'BDNF', 'system': 'BDNF expression',
        'effect': 'Depression risk',
        'CC': 'Lower expression',
        'CT': 'Intermediate',
        'TT': 'Higher BDNF expression (protective)',
    },
}

MOOD_BIPOLAR_SCHIZOPHRENIA = {
    'rs1006737': {  # CACNA1C
        'gene': 'CACNA1C', 'system': 'Calcium channel',
        'effect': 'Bipolar disorder, major depression',
        'GG': 'Typical risk',
        'AG': 'Elevated risk (~1.1x)',
        'AA': 'Higher bipolar/depression risk (~1.2x)',
    },
    'rs10994336': {  # ANK3
        'gene': 'ANK3', 'system': 'Ankyrin 3',
        'effect': 'Bipolar disorder',
        'CC': 'Typical risk',
        'CT': 'Elevated bipolar risk (1.45x)',
        'TT': 'Higher bipolar risk (2.9x)',
    },
    'rs10994397': {  # ANK3
        'gene': 'ANK3', 'system': 'Ankyrin 3',
        'effect': 'Bipolar disorder',
        'CC': 'Typical risk',
        'CT': 'Elevated risk',
        'TT': 'Higher risk',
    },
    'rs1344706': {  # ZNF804A
        'gene': 'ZNF804A', 'system': 'Zinc finger protein',
        'effect': 'Schizophrenia, bipolar',
        'GG': 'Normal risk',
        'GT': '1.1x increased schizophrenia risk',
        'TT': '1.2x increased schizophrenia risk',
        # Complement strand (some testing companies report these)
        'CC': 'Normal risk',  # Complement of GG
        'CA': '1.1x increased schizophrenia risk',  # Complement of GT
        'AC': '1.1x increased schizophrenia risk',  # Complement of TG (same as GT)
        'AA': '1.2x increased schizophrenia risk',  # Complement of TT
    },
    # rs7430407 removed - doesn't exist in SNPedia
    'rs3924999': {  # NRG1
        'gene': 'NRG1', 'system': 'Neuregulin 1',
        'effect': 'Schizophrenia association (research inconclusive, magnitude 0)',
        'GG': 'Research on schizophrenia risk is mixed/inconclusive',
        'AG': 'Research on schizophrenia risk is mixed/inconclusive',
        'AA': 'Research on schizophrenia risk is mixed/inconclusive',
    },
    'rs17512836': {  # TCF4
        'gene': 'TCF4', 'system': 'Transcription factor',
        'effect': 'Schizophrenia',
        'TT': 'Lower risk',
        'TC': 'Moderate risk (OR 1.23)',
        'CC': 'Elevated schizophrenia risk',
    },
}

MOOD_ANXIETY_STRESS = {
    'rs324420': {  # FAAH
        'gene': 'FAAH', 'system': 'Endocannabinoid breakdown',
        'effect': 'Anxiety, happiness, pain (Note: Complex effects)',
        'CC': 'Fast FAAH breakdown (typical)',
        'CA': 'Intermediate',
        'AA': 'Slow FAAH breakdown - higher pain tolerance, national happiness correlation BUT higher substance use disorder risk',
    },
    'rs53576': {  # OXTR
        'gene': 'OXTR', 'system': 'Oxytocin receptor',
        'effect': 'Social bonding, empathy, stress response',
        'GG': 'Higher sensitivity (more empathetic, handle stress well)',
        'GA': 'Intermediate',
        'AA': 'Lower sensitivity (less empathetic, higher stress)',
    },
    'rs2254298': {  # OXTR
        'gene': 'OXTR', 'system': 'Oxytocin receptor',
        'effect': 'Social behavior, autism',
        'GG': 'Typical social behavior',
        'GA': 'Moderate',
        'AA': 'Altered social behavior, autism association',
    },
    'rs237887': {  # OXTR
        'gene': 'OXTR', 'system': 'Oxytocin',
        'effect': 'Social cognition',
        'AA': 'Typical',
        'AG': 'Moderate',
        'GG': 'Altered social processing',
    },
    'rs2268498': {  # OXTR
        'gene': 'OXTR', 'system': 'Oxytocin',
        'effect': 'Social perception, facial emotion recognition',
        'CC': 'Common variant',
        'CT': 'Heterozygous',
        'TT': 'Homozygous variant (may affect social perception)',
    },
    'rs13212041': {  # HTR1B
        'gene': 'HTR1B', 'system': 'Serotonin receptor 1B',
        'effect': 'Anxiety, depression',
        'AA': 'Lower anxiety',
        'AG': 'Moderate',
        'GG': 'Higher anxiety tendency',
    },
    'rs3813034': {  # SLC6A4
        'gene': 'SLC6A4', 'system': 'Serotonin transporter',
        'effect': 'Behavior, bipolar, depression associations',
        'CC': 'Probable non-pathogenic',
        'CA': 'Heterozygous',
        'AA': 'Less common variant',
    },
    'rs1799971': {  # OPRM1
        'gene': 'OPRM1', 'system': 'Opioid receptor',
        'effect': 'Pain, addiction, alcohol response',
        'AA': 'Normal opioid response',
        'AG': 'Reduced receptor function',
        'GG': 'Lower receptor density, higher addiction risk',
    },
}

MOOD_OTHER = {
    'rs12922317': {  # SNX29
        'gene': 'SNX29', 'system': 'Sorting nexin',
        'effect': 'Possible schizophrenia association',
        'GG': 'Common genotype',
        'GA': 'Heterozygous',
        'AA': 'Less common genotype (limited research)',
    },
    'rs4675690': {  # CREB1
        'gene': 'CREB1', 'system': 'cAMP response element-binding protein',
        'effect': 'Learning, memory, schizophrenia',
        'CC': 'Typical',
        'CT': 'Moderate',
        'TT': 'Altered function',
    },
    'rs2049045': {  # BDNF
        'gene': 'BDNF', 'system': 'Brain-derived neurotrophic factor',
        'effect': 'Alzheimer risk, bipolar disorder',
        'CC': 'Common genotype',
        'CG': 'Heterozygous',
        'GG': 'Less common (limited research)',
    },
    'rs1801131_mood': {  # MTHFR A1298C
        'gene': 'MTHFR', 'system': 'Folate/neurotransmitters',
        'effect': 'Mood, neurotransmitter synthesis',
        'AA': 'Normal',
        'AC': 'Slightly reduced',
        'CC': 'Reduced function, mood impact',
    },
}


# ============================================================================
# PHYSICAL TRAITS
# ============================================================================

TRAITS_APPEARANCE = {
    'rs12913832': {  # HERC2/OCA2
        'gene': 'HERC2', 'trait': 'Eye color',
        'evidence': 'Explains ~74% of variance',
        'AA': 'Brown eyes (>95%)',
        'AG': 'Brown/green (mixed)',
        'GG': 'Blue eyes (>90%)',
    },
    'rs1800407': {  # OCA2
        'gene': 'OCA2', 'trait': 'Eye color',
        'evidence': 'Modifier',
        'GG': 'Lighter eyes (blue/gray)',
        'AG': 'Mixed (hazel/green)',
        'AA': 'Darker eyes (brown/black)',
    },
    'rs16891982': {  # SLC45A2
        'gene': 'SLC45A2', 'trait': 'Eye/skin pigmentation',
        'evidence': 'Pigmentation gene',
        'CC': 'Darker',
        'CG': 'Mixed',
        'GG': 'Lighter',
    },
    'rs1393350': {  # TYR
        'gene': 'TYR', 'trait': 'Eye color',
        'evidence': 'Tyrosinase',
        'AA': 'Lighter eyes (blue), blond hair',
        'AG': 'Lighter eyes tendency',
        'GG': 'Darker eyes (brown)',
    },
    'rs12896399': {  # SLC24A4
        'gene': 'SLC24A4', 'trait': 'Hair color',
        'evidence': 'Pigmentation',
        'GG': 'Dark hair',
        'GT': 'Mixed',
        'TT': 'Light hair',
    },
    'rs1805007': {  # MC1R
        'gene': 'MC1R', 'trait': 'Red hair',
        'evidence': 'Primary red hair gene',
        'CC': 'Non-red',
        'CT': 'Red hair carrier',
        'TT': 'Red hair likely',
    },
    'rs1805008': {  # MC1R
        'gene': 'MC1R', 'trait': 'Red hair',
        'evidence': 'Red hair variant',
        'CC': 'Non-red',
        'CT': 'Carrier',
        'TT': 'Red hair',
    },
    'rs1805009': {  # MC1R
        'gene': 'MC1R', 'trait': 'Red hair',
        'evidence': 'Red hair variant',
        'CC': 'Non-red',
        'CT': 'Carrier',
        'TT': 'Red hair',
    },
    'rs11547464': {  # MC1R
        'gene': 'MC1R', 'trait': 'Red hair',
        'evidence': 'Red hair variant',
        'AA': 'Red hair likely',
        'AG': 'Red hair carrier',
        'GG': 'Non-red hair',
    },
    'rs1110400': {  # MC1R
        'gene': 'MC1R', 'trait': 'Melanoma risk',
        'evidence': 'MC1R region - ClinVar pathogenic',
        'CC': 'Higher melanoma risk',
        'CT': 'Moderate risk',
        'TT': 'Lower risk',
    },
    'rs17646946': {  # TCHH
        'gene': 'TCHH', 'trait': 'Hair texture',
        'evidence': 'Trichohyalin',
        'AA': 'Straight hair',
        'AG': 'Mixed',
        'GG': 'Curly/wavy',
    },
    'rs3827760': {  # EDAR
        'gene': 'EDAR', 'trait': 'Hair texture/thickness',
        'evidence': 'East Asian variant',
        'CC': 'Straight, thick hair',
        'CT': 'Mixed',
        'TT': 'Variable',
    },
    'rs1426654': {  # SLC24A5
        'gene': 'SLC24A5', 'trait': 'Skin color',
        'evidence': 'Major determinant',
        'AA': 'Lighter skin (European ancestry)',
        'AG': 'Intermediate',
        'GG': 'Darker skin (African/Asian ancestry)',
    },
    'rs12203592': {  # IRF4
        'gene': 'IRF4', 'trait': 'Freckles',
        'evidence': 'Freckling tendency',
        'CC': 'No freckles',
        'CT': 'Some freckles',
        'TT': 'More freckles',
    },
}

TRAITS_BODY = {
    'rs1815739': {  # ACTN3
        'gene': 'ACTN3', 'trait': 'Muscle fiber type',
        'evidence': 'Olympic athlete studies',
        'CC': 'Fast-twitch (RR - sprinter)',
        'CT': 'Mixed (RX)',
        'TT': 'Endurance (XX - no alpha-actinin-3)',
    },
    'rs8111989': {  # NFIA
        'gene': 'NFIA', 'trait': 'Endurance',
        'evidence': 'Endurance performance',
        'AA': 'Endurance advantage',
        'AG': 'Mixed',
        'GG': 'Typical',
    },
    'rs1042713': {  # ADRB2
        'gene': 'ADRB2', 'trait': 'Asthma/inhaler response',
        'evidence': 'Beta-2 adrenergic receptor',
        'AA': 'Arg/Arg - 1.7x risk that pediatric inhaler use may worsen asthma',
        'AG': 'Gly/Arg - 1.3x risk that pediatric inhaler use may worsen asthma',
        'GG': 'Gly/Gly - normal (reference)',
    },
    'rs1042714': {  # ADRB2
        'gene': 'ADRB2', 'trait': 'Multiple health associations (Q27E)',
        'evidence': 'Beta-2 adrenergic receptor',
        'CC': 'Gln/Gln - normal/baseline risk',
        'CG': 'Gln/Glu - complex (increased thromboembolism, autism risk)',
        'GG': 'Glu/Glu - complex (higher thromboembolism, autism, stroke; lower diabetes risk)',
    },
    'rs6570507': {  # ADGRG6
        'gene': 'ADGRG6', 'trait': 'Height',
        'evidence': 'Growth factor',
        'CC': 'Shorter tendency',
        'CT': 'Average',
        'TT': 'Taller tendency',
    },
    'rs1042725': {  # HMGA2
        'gene': 'HMGA2', 'trait': 'Height',
        'evidence': '~0.4cm per C allele',
        'TT': 'Shorter',
        'TC': 'Average (+0.4cm)',
        'CC': 'Taller (+0.8cm)',
    },
    'rs17822931': {  # ABCC11
        'gene': 'ABCC11', 'trait': 'Earwax type/body odor',
        'evidence': 'Well-established',
        'CC': 'Wet earwax, more body odor',
        'CT': 'Wet (carrier for dry)',
        'TT': 'Dry earwax, less body odor',
    },
}

TRAITS_METABOLISM = {
    'rs4988235': {  # LCT
        'gene': 'LCT', 'trait': 'Lactose tolerance',
        'evidence': 'Primary European variant',
        'CC': 'Lactose intolerant (low lactase activity)',
        'CT': 'Likely tolerant (intermediate activity)',
        'TT': 'Lactose tolerant (can digest milk)',
    },
    'rs182549': {  # LCT
        'gene': 'LCT', 'trait': 'Lactose intolerance',
        'evidence': 'Additional variant',
        'CC': 'Intolerant',
        'CT': 'Likely tolerant',
        'TT': 'Tolerant',
    },
    'rs713598': {  # TAS2R38
        'gene': 'TAS2R38', 'trait': 'Bitter taste (PTC)',
        'evidence': 'Taste receptor',
        'CC': 'Taster',
        'CG': 'Taster',
        'GG': 'Non-taster',
    },
    'rs1726866': {  # TAS2R38
        'gene': 'TAS2R38', 'trait': 'Bitter taste',
        'evidence': 'Taste',
        'CC': 'Taster',
        'CT': 'Taster',
        'TT': 'Non-taster',
    },
    'rs671': {  # ALDH2
        'gene': 'ALDH2', 'trait': 'Alcohol flush',
        'evidence': 'Asian flush',
        'GG': 'Normal metabolism',
        'GA': 'Reduced (flush reaction)',
        'AA': 'Very reduced (strong flush)',
    },
    'rs1229984': {  # ADH1B
        'gene': 'ADH1B', 'trait': 'Alcohol metabolism',
        'evidence': 'Alcohol dehydrogenase',
        'AA': 'Typical',
        'AG': 'Faster',
        'GG': 'Fast (protective against alcoholism)',
    },
    'rs762551': {  # CYP1A2
        'gene': 'CYP1A2', 'trait': 'Caffeine metabolism',
        'evidence': 'Main metabolizer',
        'AA': 'Fast metabolizer (especially in smokers/heavy coffee drinkers)',
        'AC': 'Normal metabolism',
        'CC': 'Normal metabolism',
    },
    'rs2472297': {  # Between CYP1A1 and CYP1A2
        'gene': 'Between CYP1A1/CYP1A2', 'trait': 'Caffeine',
        'evidence': 'Metabolism (intergenic region)',
        'CC': 'Fast',
        'CT': 'Intermediate',
        'TT': 'Slower',
    },
    'rs16969968': {  # CHRNA5
        'gene': 'CHRNA5', 'trait': 'Nicotine dependence',
        'evidence': 'Smoking behavior',
        'AA': 'Higher dependence risk',
        'AG': 'Moderate',
        'GG': 'Lower risk',
    },
    'rs1051730': {  # CHRNA3
        'gene': 'CHRNA3', 'trait': 'Nicotine dependence',
        'evidence': 'Smoking, lung cancer',
        'CC': 'Lower risk (fewer cigarettes if smoker)',
        'CT': 'Moderate (1.3x lung cancer risk)',
        'TT': 'Higher risk (1.8x lung cancer risk)',
    },
    'rs72921001': {  # OR6A2
        'gene': 'OR6A2', 'trait': 'Cilantro taste',
        'evidence': 'Olfactory receptor',
        'AA': 'Normal taste',
        'AG': 'May taste soapy',
        'GG': 'Soapy taste',
    },
    'rs4481887': {  # Intergenic
        'gene': 'Intergenic', 'trait': 'Asparagus odor detection',
        'evidence': 'Olfactory',
        'AA': 'Can smell',
        'AG': 'Can smell',
        'GG': 'Cannot smell',
    },
}

TRAITS_SLEEP_CIRCADIAN = {
    'rs1801260': {  # CLOCK
        'gene': 'CLOCK', 'trait': 'Circadian rhythm',
        'evidence': 'Clock gene',
        'CC': 'Strong evening preference (delayed sleep ~79min, less sleep ~75min)',
        'CT': 'Moderate evening preference',
        'TT': 'Normal/typical sleep pattern',
    },
    'rs73598374': {  # ADA
        'gene': 'ADA', 'trait': 'Deep sleep',
        'evidence': 'Sleep quality',
        'AA': 'Less deep sleep',
        'AG': 'Normal',
        'GG': 'More deep sleep',
    },
    'rs2292912': {  # CRY2
        'gene': 'CRY2', 'trait': 'Sleep timing',
        'effect': 'Morning/evening chronotype',
        'CC': 'Morning person',
        'CT': 'Intermediate',
        'TT': 'Evening person',
    },
}

TRAITS_VITAMINS = {
    'rs4588': {  # GC
        'gene': 'GC', 'trait': 'Vitamin D levels',
        'evidence': 'Vitamin D binding protein - affects 25(OH)D levels',
        'CC': 'Higher vitamin D levels (C allele associated with ~10-15% higher 25(OH)D)',
        'CA': 'Intermediate vitamin D levels',
        'AA': 'Lower vitamin D levels (A allele associated with ~10-15% lower 25(OH)D)',
    },
    'rs2282679': {  # GC
        'gene': 'GC', 'trait': 'Vitamin D',
        'evidence': 'Vitamin D',
        'GG': 'Lower levels',
        'GT': 'Moderate',
        'TT': 'Higher levels',
    },
    'rs601338': {  # FUT2
        'gene': 'FUT2', 'trait': 'B12 levels',
        'evidence': 'Secretor status',
        'AA': 'Non-secretor (higher plasma B12 levels)',
        'AG': 'Intermediate',
        'GG': 'Secretor (lower plasma B12 levels)',
    },
    'rs4654748': {  # NBPF3
        'gene': 'NBPF3', 'trait': 'Vitamin B6',
        'evidence': 'B6 blood concentration',
        'CC': 'Lower B6 levels (-1.45 to -2.90 ng/mL)',
        'CT': 'Moderate B6 levels',
        'TT': 'Higher B6 levels',
    },
}

TRAITS_ATHLETIC = {
    'rs1815739_athletic': {  # ACTN3
        'gene': 'ACTN3', 'trait': 'Athletic performance',
        'evidence': 'R577X - Olympic athlete studies',
        'CC': 'RR - Power/sprint athlete advantage',
        'CT': 'RX - Mixed performance',
        'TT': 'XX - Endurance athlete advantage',
    },
    'rs8192678': {  # PPARGC1A
        'gene': 'PPARGC1A', 'trait': 'Aerobic capacity, V̇O2max',
        'evidence': 'PGC-1α - mitochondrial biogenesis',
        'CC': 'Better aerobic capacity',
        'CT': 'Intermediate',
        'TT': 'Lower V̇O2max response to training',
    },
    'rs1572312': {  # NFIA-AS2
        'gene': 'NFIA-AS2', 'trait': 'Sports genomics marker',
        'evidence': 'Limited research (alleles: A/C)',
        'CC': 'Common genotype',
        'CA': 'Heterozygous',
        'AA': 'Alternative genotype',
    },
    'rs699': {  # AGT
        'gene': 'AGT', 'trait': 'Blood pressure response to exercise',
        'evidence': 'Angiotensinogen',
        'TT': 'Lower blood pressure',
        'CT': 'Intermediate',
        'CC': 'Higher blood pressure, power athlete association',
    },
}

TRAITS_SKIN_AGING = {
    'rs1800975': {  # MMP1
        'gene': 'MMP1', 'trait': 'Skin aging, wrinkles',
        'evidence': 'Collagen breakdown',
        'GG': 'Lower collagen breakdown',
        'GA': 'Moderate',
        'AA': 'Higher collagen breakdown, premature aging',
    },
    'rs1800012': {  # COL1A1
        'gene': 'COL1A1', 'trait': 'Collagen strength, bone density',
        'evidence': 'Type I collagen Sp1 polymorphism (meta-analysis: PMC5650439)',
        'GG': 'Normal/optimal collagen production, highest bone density, lowest fracture risk',
        'GT': 'Intermediate - moderately reduced bone density',
        'TT': 'Increased fracture risk (OR=1.78-2.97), lower bone density, imbalanced collagen ratio (α1:α2 = 2.3:1 vs normal 2:1)',
    },
    'rs1126809': {  # TYR
        'gene': 'TYR', 'trait': 'Skin pigmentation, melanoma risk',
        'evidence': 'Tyrosinase (R402Q)',
        'GG': 'Arg/Arg - normal melanoma risk',
        'GA': 'Arg/Gln - slight increase in skin cancer risk',
        'AA': 'Gln/Gln - slight increase in skin cancer risk',
    },
}

TRAITS_TASTE_EXTENDED = {
    'rs307355': {  # TAS1R3
        'gene': 'TAS1R3', 'trait': 'Sweet taste perception',
        'evidence': 'Sweet receptor',
        'CC': 'Normal sweet perception',
        'CT': 'Reduced sweet sensitivity',
        'TT': 'Lower sweet taste sensitivity',
    },
    'rs35874116': {  # TAS1R2
        'gene': 'TAS1R2', 'trait': 'Sweet taste',
        'evidence': 'Sweet receptor',
        'AA': 'Normal',
        'AG': 'Altered sweet perception',
        'GG': 'Reduced sweet taste',
    },
    'rs846664': {  # TRPV1
        'gene': 'TRPV1', 'trait': 'Spicy food tolerance',
        'evidence': 'Capsaicin receptor',
        'AA': 'Higher spice tolerance',
        'AG': 'Moderate',
        'GG': 'Lower spice tolerance',
    },
}

TRAITS_THERMOGENESIS = {
    'rs1801282_thermo': {  # PPARG
        'gene': 'PPARG', 'trait': 'Fat storage, thermogenesis',
        'evidence': 'Pro12Ala affects energy expenditure',
        'CC': 'Higher fat storage tendency',
        'CG': 'Moderate',
        'GG': 'Better thermogenesis, lower fat storage',
    },
    'rs1800592': {  # UCP1
        'gene': 'UCP1', 'trait': 'Thermogenesis, cold tolerance',
        'evidence': 'Uncoupling protein - brown fat',
        'GG': 'Lower thermogenesis',
        'GA': 'Moderate',
        'AA': 'Higher thermogenesis, better cold tolerance',
    },
    'rs659366': {  # UCP2
        'gene': 'UCP2', 'trait': 'Energy expenditure, obesity',
        'evidence': 'Mitochondrial uncoupling',
        'AA': 'Higher energy expenditure',
        'AC': 'Moderate',
        'CC': 'Lower energy expenditure, obesity risk',
    },
}


# ============================================================================
# ANALYSIS ENGINE
# ============================================================================

def get_genotype_interpretation(genotype, data):
    """Get interpretation handling reversed order and strand orientation"""
    # Direct match
    if genotype in data:
        return data[genotype]

    # Try reversed (TC -> CT)
    reversed_gt = genotype[::-1]
    if reversed_gt in data:
        return data[reversed_gt]

    # Try complement strand
    complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    complement = ''.join([complement_map.get(a, a) for a in genotype])
    if complement in data:
        return data[complement]

    # Try reversed complement
    reversed_complement = complement[::-1]
    if reversed_complement in data:
        return data[reversed_complement]

    return f"Unknown genotype: {genotype}"


def analyze_category(snps_obj, markers, category_name):
    """Generic analyzer for any marker category"""
    results = []
    for rsid, data in markers.items():
        snp_data = snps_obj.snps[snps_obj.snps.index == rsid]
        if len(snp_data) > 0 and not snp_data['genotype'].isna().all():
            genotype = str(snp_data['genotype'].iloc[0])
            gene = data.get('gene', '')
            interpretation = get_genotype_interpretation(genotype, data)

            result = {
                'rsid': rsid,
                'gene': gene,
                'genotype': genotype,
                'interpretation': interpretation,
            }

            # Add category-specific fields
            for key in ['drug', 'condition', 'trait', 'system', 'effect', 'evidence']:
                if key in data:
                    result[key] = data[key]

            results.append(result)

    return results


def get_interpretation_color(interpretation):
    """Determine color based on interpretation text"""
    interpretation_lower = interpretation.lower()

    # Good/protective indicators - CHECK FIRST (highest priority)
    good_keywords = [
        'protective', 'better', 'fast metabol',
        'advantage', 'more empathetic',
        'resilient', 'warrior', 'tolerant', 'optimal', 'enhanced',
        'improved', 'better memory maintenance'
    ]

    # Protective risk phrases (contain "risk" but are good)
    protective_risk_phrases = [
        'lower risk', 'reduced risk', 'decreased risk',
        'lower alzheimer', 'reduced alzheimer', 'decreased alzheimer',
        'lower ad risk', 'reduced ad risk',
        'typical risk', 'normal risk'
    ]

    # Check for protective/good phrases first
    if any(keyword in interpretation_lower for keyword in good_keywords):
        return Colors.GREEN

    if any(phrase in interpretation_lower for phrase in protective_risk_phrases):
        return Colors.GREEN

    # Bad/risk indicators - RED
    bad_keywords = [
        'elevated risk', 'higher risk', 'high risk', 'increased risk',
        'poor metabol', 'reduced', 'deficiency', 'carrier',
        'slow metabol', 'worse', 'poorer', 'depression risk', 'anxiety risk',
        'addiction risk', 'disease risk', 'overload', 'mutation',
        'homozygous', 'much higher', 'strong flush', 'less empathetic',
        'reduced receptor', 'lower d2', 'lower expression', 'altered',
        'reduced function', 'impaired', 'susceptible', 'faster decline',
        'accelerated cognitive decline', 'poorer memory'
    ]

    # Neutral/moderate indicators - YELLOW
    moderate_keywords = [
        'intermediate', 'moderate', 'mixed', 'carrier', 'typical',
        'normal', 'heterozygous', 'may need', 'slightly'
    ]

    # Check for bad
    if any(keyword in interpretation_lower for keyword in bad_keywords):
        return Colors.RED

    # Then check for moderate
    if any(keyword in interpretation_lower for keyword in moderate_keywords):
        return Colors.YELLOW

    # Default to white/no color
    return ''


def print_results(results, title, show_evidence=True):
    """Print results in a clean format"""
    if not results:
        return

    print(f"\n{Colors.BOLD}{Colors.YELLOW}{title}:{Colors.END}")
    print(f"{Colors.GRAY}{'-' * 80}{Colors.END}")

    for r in results:
        gene_display = f"{Colors.BOLD}{Colors.GREEN}{r['gene']}{Colors.END}" if r['gene'] else ""
        rsid_display = f"{Colors.GRAY}({r['rsid']}){Colors.END}"

        # Build the header line
        header_parts = [gene_display, rsid_display]
        if 'drug' in r:
            header_parts.insert(1, f"{Colors.ORANGE}[{r['drug']}]{Colors.END}")
        elif 'condition' in r:
            header_parts.insert(1, f"{Colors.ORANGE}[{r['condition']}]{Colors.END}")
        elif 'trait' in r:
            header_parts.insert(1, f"{Colors.BLUE}[{r['trait']}]{Colors.END}")
        elif 'system' in r:
            header_parts.insert(1, f"{Colors.HEADER}[{r['system']}]{Colors.END}")

        print(f"\n  {' '.join(header_parts)}: {Colors.BOLD}{r['genotype']}{Colors.END}")

        if 'effect' in r:
            print(f"  {Colors.GRAY}Effect:{Colors.END} {r['effect']}")

        if show_evidence and 'evidence' in r:
            print(f"  {Colors.GRAY}Evidence:{Colors.END} {r['evidence']}")

        # Color-code the interpretation based on content
        interp_color = get_interpretation_color(r['interpretation'])
        print(f"  {Colors.BOLD}→{Colors.END} {interp_color}{r['interpretation']}{Colors.END}")
        print(f"  {Colors.BLUE}🔗 https://www.snpedia.com/index.php/{r['rsid']}{Colors.END}")


def analyze_dna(filepath):
    """Main analysis function"""

    print(f"\n{Colors.BOLD}{Colors.CYAN}{'='*80}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN} COMPREHENSIVE DNA ANALYSIS{Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}{'='*80}{Colors.END}")
    print(f"\n{Colors.YELLOW}Loading DNA data...{Colors.END}\n")

    s = snps.SNPs(filepath)

    # Basic Info
    print_section("BASIC INFORMATION")
    print(f"{Colors.BOLD}Source:{Colors.END} {Colors.GREEN}{s.source}{Colors.END}")
    print(f"{Colors.BOLD}Total SNPs:{Colors.END} {Colors.GREEN}{s.count:,}{Colors.END}")
    print(f"{Colors.BOLD}Build:{Colors.END} {Colors.GREEN}{s.build}{Colors.END}")
    print(f"{Colors.BOLD}Sex:{Colors.END} {Colors.GREEN}{s.sex}{Colors.END}")

    # Tier 1: Clinical Grade
    print_section("TIER 1: CLINICAL GRADE MARKERS")
    print(f"{Colors.GRAY}FDA recognized | ACMG reportable | Direct clinical actionability{Colors.END}\n")

    results = analyze_category(s, TIER1_PHARMACOGENOMICS, "Pharmacogenomics")
    print_results(results, "PHARMACOGENOMICS (FDA Recognized)", show_evidence=True)

    results = analyze_category(s, TIER1_CARDIOVASCULAR, "Cardiovascular")
    print_results(results, "CARDIOVASCULAR & IRON (ACMG Reportable)", show_evidence=True)

    # Tier 2: Health Markers
    print_section("TIER 2: HEALTH & DISEASE RISK")
    print(f"{Colors.GRAY}GWAS validated | Strong research evidence{Colors.END}\n")

    for category_dict, category_name in [
        (HEALTH_METABOLIC, "Metabolic & Diabetes"),
        (HEALTH_CARDIOVASCULAR_EXTENDED, "Cardiovascular Extended"),
        (HEALTH_CANCER, "Cancer Risk Markers"),
        (HEALTH_COVID19, "COVID-19 & Blood Type"),
        (HEALTH_AUTOIMMUNE, "Autoimmune & Inflammation"),
        (HEALTH_BONE_KIDNEY, "Bone, Kidney & Vascular"),
        (LONGEVITY_COGNITIVE_AGING, "Longevity & Cognitive Aging ⭐"),
        (HEALTH_COGNITIVE, "Cognitive Function & Memory"),
        (HEALTH_PAIN, "Pain Sensitivity"),
        (HEALTH_ADDICTION, "Addiction Susceptibility"),
    ]:
        results = analyze_category(s, category_dict, category_name)
        print_results(results, category_name.upper(), show_evidence=False)

        # Print personalized APOE + FOXO3 interpretation after longevity section
        if category_dict is LONGEVITY_COGNITIVE_AGING:
            interpret_apoe_foxo3(s)

    # Mood & Mental Health
    print_section("MOOD, MENTAL HEALTH & NEUROTRANSMITTERS")

    for category_dict, category_name in [
        (MOOD_SEROTONIN, "Serotonin System"),
        (MOOD_DOPAMINE, "Dopamine System"),
        (MOOD_NEUROPLASTICITY, "Neuroplasticity (BDNF)"),
        (MOOD_BIPOLAR_SCHIZOPHRENIA, "Bipolar & Schizophrenia"),
        (MOOD_ANXIETY_STRESS, "Anxiety, Stress & Social Bonding"),
        (MOOD_OTHER, "Other Neurotransmitter Systems"),
    ]:
        results = analyze_category(s, category_dict, category_name)
        print_results(results, category_name.upper(), show_evidence=False)

    # Physical Traits
    print_section("PHYSICAL TRAITS & CHARACTERISTICS")

    for category_dict, category_name in [
        (TRAITS_APPEARANCE, "Appearance (Eyes, Hair, Skin)"),
        (TRAITS_BODY, "Body & Performance"),
        (TRAITS_ATHLETIC, "Athletic Performance"),
        (TRAITS_SKIN_AGING, "Skin Aging & Collagen"),
        (TRAITS_METABOLISM, "Taste, Smell & Substance Metabolism"),
        (TRAITS_TASTE_EXTENDED, "Taste Perception (Sweet, Spicy)"),
        (TRAITS_SLEEP_CIRCADIAN, "Sleep & Circadian Rhythm"),
        (TRAITS_VITAMINS, "Vitamin Metabolism"),
        (TRAITS_THERMOGENESIS, "Thermogenesis & Energy Expenditure"),
    ]:
        results = analyze_category(s, category_dict, category_name)
        print_results(results, category_name.upper(), show_evidence=False)

    # Check CLI arguments if they exist
    cli_args = globals().get('CLI_ARGS', None)

    # Athletic Polygenic Score (prints its own header)
    if not (cli_args and cli_args.no_athletic):
        analyze_athletic_performance(s)

    # Mitochondrial DNA section removed - not useful without haplogroup analysis

    # Summary
    print_section("ANALYSIS COMPLETE")
    print(f"{Colors.GREEN}✓ All processing done locally - no data transmitted{Colors.END}")

    # Count markers found
    total_markers = 0
    for marker_dict in [
        TIER1_PHARMACOGENOMICS, TIER1_CARDIOVASCULAR,
        HEALTH_METABOLIC, HEALTH_CARDIOVASCULAR_EXTENDED, HEALTH_CANCER,
        HEALTH_COVID19, HEALTH_AUTOIMMUNE, HEALTH_BONE_KIDNEY, HEALTH_LONGEVITY,
        HEALTH_COGNITIVE, HEALTH_PAIN, HEALTH_ADDICTION,
        MOOD_SEROTONIN, MOOD_DOPAMINE, MOOD_NEUROPLASTICITY,
        MOOD_BIPOLAR_SCHIZOPHRENIA, MOOD_ANXIETY_STRESS, MOOD_OTHER,
        TRAITS_APPEARANCE, TRAITS_BODY, TRAITS_ATHLETIC, TRAITS_SKIN_AGING,
        TRAITS_METABOLISM, TRAITS_TASTE_EXTENDED, TRAITS_SLEEP_CIRCADIAN,
        TRAITS_VITAMINS, TRAITS_THERMOGENESIS,
    ]:
        results = analyze_category(s, marker_dict, "")
        total_markers += len(results)

    print(f"\n{Colors.BOLD}Total markers analyzed:{Colors.END} {Colors.GREEN}{total_markers}{Colors.END}")
    print(f"\n{Colors.BOLD}{Colors.YELLOW}Resources:{Colors.END}")
    print(f"  {Colors.BLUE}• SNPedia.com - Detailed SNP information{Colors.END}")
    print(f"  {Colors.BLUE}• PharmGKB.org - Pharmacogenomics{Colors.END}")
    print(f"  {Colors.BLUE}• ClinVar - Clinical variants{Colors.END}")
    print(f"  {Colors.BLUE}• CPIC guidelines: cpicpgx.org{Colors.END}")

    return s


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Analyze AncestryDNA raw data for health, traits, and ancestry insights',
        epilog='Examples:\n'
               '  %(prog)s                           # Use default AncestryDNA.txt\n'
               '  %(prog)s mydata.txt                # Analyze specific file\n'
               '  %(prog)s data.txt --no-cognitive   # Skip cognitive score\n'
               '  %(prog)s data.txt --no-athletic    # Skip athletic score\n'
               '  %(prog)s data.txt --health-only    # Only health analysis\n'
               '  %(prog)s data.txt --quick          # Skip polygenic scores\n',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('filepath', nargs='?', default='AncestryDNA.txt',
                        help='Path to AncestryDNA raw data file (default: AncestryDNA.txt)')
    parser.add_argument('--no-athletic', action='store_true',
                        help='Skip athletic performance analysis')
    parser.add_argument('--health-only', action='store_true',
                        help='Only run health-related analyses')
    parser.add_argument('--quick', action='store_true',
                        help='Quick mode: skip polygenic scores (athletic)')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Show detailed output and statistics')

    args = parser.parse_args()

    # Handle quick mode and health-only mode
    if args.quick:
        args.no_athletic = True
    if args.health_only:
        args.no_athletic = True

    try:
        # Store args globally so analyze_dna function can access them
        globals()['CLI_ARGS'] = args

        snps_data = analyze_dna(args.filepath)
        print(f"\n{Colors.BOLD}{Colors.GREEN}{'='*80}{Colors.END}")
        print(f"{Colors.BOLD}{Colors.GREEN} DONE{Colors.END}")
        print(f"{Colors.BOLD}{Colors.GREEN}{'='*80}{Colors.END}\n")

    except FileNotFoundError:
        print(f"{Colors.RED}Error: Could not find file '{args.filepath}'{Colors.END}")
        print(f"{Colors.YELLOW}Try: python {parser.prog} --help{Colors.END}")
        sys.exit(1)
    except Exception as e:
        print(f"{Colors.RED}Error: {e}{Colors.END}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)
