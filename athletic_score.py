#!/usr/bin/env python3
"""
Training Bias & Durability Profile (Limited SNP Panel)

Analyzes genetic variants across 4 orthogonal dimensions:
1. OUTPUT INDEX - Force/speed/power expression (what wins races/lifts)
2. DURABILITY INDEX - Work capacity, injury resistance (training tolerance)
3. ADAPTATION INDEX - Skill acquisition, neuroplasticity (learning rate)
4. METABOLIC INDEX - Body composition, energy efficiency (health-adjacent, indirect)

NOT a general "athletic performance" predictor.
Each index measures different causal mechanisms - aggregation destroys meaning.
"""

import snps
from typing import Dict, Tuple, Optional, List
from dataclasses import dataclass
import math


@dataclass
class AthleticSNP:
    """Athletic performance SNP with mechanistic categorization"""
    rsid: str
    gene: str
    name: str
    index: str  # 'output', 'durability', 'adaptation', 'metabolic'
    favorable_alleles: List[str]  # Alleles that improve the measured trait
    unfavorable_alleles: List[str]  # Alleles that impair the measured trait
    effect_description: str
    evidence: str
    weight: float  # Brutal reweighting: ACE=1.0, ACTN3=0.5, rest 0.05-0.4
    notes: Optional[str] = None
    snpedia_url: Optional[str] = None


# Athletic Performance SNPs Database
# Brutally reweighted - only what matters gets weight
ATHLETIC_SNPS = {
    # ============================================================================
    # OUTPUT INDEX - Force/speed/power expression
    # ============================================================================

    'rs4343': AthleticSNP(
        rsid='rs4343', gene='ACE', name='ACE I/D proxy',
        index='output',
        favorable_alleles=['G'],  # D allele (power)
        unfavorable_alleles=['A'],  # I allele
        effect_description='ACE D allele: force/power expression',
        evidence='Elite athlete studies, replicated in power sports',
        weight=1.0,  # THE big one
        notes='G=D (deletion), A=I (insertion). D dominates power athletes.',
        snpedia_url='https://www.snpedia.com/index.php/Rs4343'
    ),

    'rs1815739': AthleticSNP(
        rsid='rs1815739', gene='ACTN3', name='R577X',
        index='output',
        favorable_alleles=['C'],  # R allele (fast-twitch)
        unfavorable_alleles=['T'],  # X allele (no alpha-actinin-3)
        effect_description='Alpha-actinin-3 in fast-twitch fibers',
        evidence='Olympic athlete studies, meta-analyses',
        weight=0.5,  # Second tier
        notes='CC (RR) = fast-twitch; TT (XX) = no alpha-actinin-3',
        snpedia_url='https://www.snpedia.com/index.php/Rs1815739'
    ),

    'rs17602729': AthleticSNP(
        rsid='rs17602729', gene='AMPD1', name='Q12X',
        index='output',
        favorable_alleles=['G'],  # Normal enzyme (Q12)
        unfavorable_alleles=['A'],  # Deficiency (stop codon)
        effect_description='AMP deaminase - muscle energy burst',
        evidence='OR=2.17 for power deficit with deficiency (PMC12152022)',
        weight=0.4,
        notes='GG = normal enzyme; AA = deficiency, power drop',
        snpedia_url='https://www.snpedia.com/index.php/Rs17602729'
    ),

    'rs2306862': AthleticSNP(
        rsid='rs2306862', gene='LRP5', name='LRP5 Wnt',
        index='output',
        favorable_alleles=['C'],  # Lean mass, strength
        unfavorable_alleles=['T'],
        effect_description='LRP5 Wnt signaling, lean mass, bone density',
        evidence='Strong lean-mass GWAS signal, grip strength',
        weight=0.15,
        notes='CC = higher lean mass, better grip strength',
        snpedia_url='https://www.snpedia.com/index.php/Rs2306862'
    ),

    'rs2228570': AthleticSNP(
        rsid='rs2228570', gene='VDR', name='FokI',
        index='output',
        favorable_alleles=['A'],  # F allele (more active VDR)
        unfavorable_alleles=['G'],  # f allele
        effect_description='Vitamin D receptor - muscle fiber composition',
        evidence='Replicated links to strength, jump height, injury risk',
        weight=0.3,
        notes='AA (FF) = more active VDR, type-II fiber advantage',
        snpedia_url='https://www.snpedia.com/index.php/Rs2228570'
    ),

    'rs1805086': AthleticSNP(
        rsid='rs1805086', gene='MSTN', name='K153R',
        index='output',
        favorable_alleles=['G'],  # R allele (reduced myostatin)
        unfavorable_alleles=['A'],  # Normal myostatin
        effect_description='Myostatin reduction, muscle mass regulation',
        evidence='R allele OR=2.02 in strength athletes (PMC3024427)',
        weight=0.3,
        notes='GG (RR) = reduced myostatin, more muscle mass (ultra-rare ~3-4%)',
        snpedia_url='https://www.snpedia.com/index.php/Rs1805086'
    ),

    # ============================================================================
    # DURABILITY INDEX - Work capacity, injury resistance
    # ============================================================================

    'rs12722': AthleticSNP(
        rsid='rs12722', gene='COL5A1', name='Type V collagen',
        index='durability',
        favorable_alleles=['T'],  # Optimal tendon stiffness
        unfavorable_alleles=['C'],
        effect_description='Collagen V structure, tendon stiffness',
        evidence='Running economy, injury risk studies',
        weight=0.4,
        notes='TT = optimal tendon stiffness, lower injury risk',
        snpedia_url='https://www.snpedia.com/index.php/Rs12722'
    ),

    'rs1800012': AthleticSNP(
        rsid='rs1800012', gene='COL1A1', name='Type I collagen Sp1',
        index='durability',
        favorable_alleles=['G'],  # Normal collagen
        unfavorable_alleles=['T'],  # Increased fracture risk
        effect_description='Collagen I production, bone mineral density',
        evidence='T allele OR=1.78 for fracture, lower BMD (meta-analysis)',
        weight=0.4,
        notes='GG = normal; TT = increased fracture risk, lower bone density',
        snpedia_url='https://www.snpedia.com/index.php/Rs1800012'
    ),

    'rs42524': AthleticSNP(
        rsid='rs42524', gene='COL1A2', name='Type I collagen α2',
        index='durability',
        favorable_alleles=['G'],  # Protective
        unfavorable_alleles=['A'],  # Injury risk
        effect_description='Type I collagen α2 chain, complements COL1A1',
        evidence='Injury risk, durability associations',
        weight=0.4,
        notes='Complements COL1A1 for connective tissue integrity'
    ),

    'rs970547': AthleticSNP(
        rsid='rs970547', gene='COL12A1', name='Type XII collagen',
        index='durability',
        favorable_alleles=['A'],  # Lower injury risk
        unfavorable_alleles=['G'],  # Higher risk
        effect_description='Tendon stiffness, ACL/soft tissue injury',
        evidence='Soft tissue injury associations, tendon mechanics',
        weight=0.4,
        notes='AA = better tendon mechanics, lower soft tissue injury risk'
    ),

    'rs4880': AthleticSNP(
        rsid='rs4880', gene='SOD2', name='Ala16Val',
        index='durability',
        favorable_alleles=['C'],  # Ala (better import, higher activity)
        unfavorable_alleles=['T'],  # Val (reduced activity)
        effect_description='Mitochondrial superoxide dismutase, ROS protection',
        evidence='Ala 30-40% more active than Val, mitochondrial import',
        weight=0.15,
        notes='CC (Ala/Ala) = better mitochondrial protection, recovery',
        snpedia_url='https://www.snpedia.com/index.php/Rs4880'
    ),

    'rs6721961': AthleticSNP(
        rsid='rs6721961', gene='NFE2L2', name='NRF2 -617C>A',
        index='durability',
        favorable_alleles=['C'],  # Normal Nrf2
        unfavorable_alleles=['A'],  # Reduced Nrf2
        effect_description='NRF2 antioxidant response regulator',
        evidence='Recovery tolerance, fatigue resistance, oxidative stress',
        weight=0.15,
        notes='CC = better antioxidant response, faster recovery',
        snpedia_url='https://www.snpedia.com/index.php/Rs6721961'
    ),

    'rs1800255': AthleticSNP(
        rsid='rs1800255', gene='COL3A1', name='Ala698Thr',
        index='durability',
        favorable_alleles=['G'],  # Normal collagen structure
        unfavorable_alleles=['A'],  # Disrupted triple helix
        effect_description='Type III collagen structure',
        evidence='A allele disrupts triple helix, OR=4.79 for injury',
        weight=0.3,
        notes='GG = normal structure; AA = disrupted, higher injury risk',
        snpedia_url='https://www.snpedia.com/index.php/Rs1800255'
    ),

    'rs2228145': AthleticSNP(
        rsid='rs2228145', gene='IL6R', name='Asp358Ala',
        index='durability',
        favorable_alleles=['C'],  # 358Ala (anti-inflammatory)
        unfavorable_alleles=['A'],  # Normal IL-6R
        effect_description='IL-6 receptor shedding, inflammation',
        evidence='C allele lowers CRP, anti-inflammatory',
        weight=0.1,
        notes='CC = lower inflammation, better recovery',
        snpedia_url='https://www.snpedia.com/index.php/Rs2228145'
    ),

    'rs1205': AthleticSNP(
        rsid='rs1205', gene='CRP', name='C-reactive protein',
        index='durability',
        favorable_alleles=['T'],  # Lowers CRP 20%
        unfavorable_alleles=['C'],  # Higher CRP
        effect_description='Baseline inflammation levels',
        evidence='Each T allele lowers CRP by 20%',
        weight=0.05,
        notes='TT = lowest inflammation, faster recovery',
        snpedia_url='https://www.snpedia.com/index.php/Rs1205'
    ),

    'rs4253778': AthleticSNP(
        rsid='rs4253778', gene='PPARA', name='PPAR-alpha',
        index='durability',
        favorable_alleles=['G'],  # Better fat oxidation
        unfavorable_alleles=['C'],
        effect_description='Fat oxidation capacity, work tolerance',
        evidence='Endurance athlete enrichment, metabolic capacity',
        weight=0.1,
        notes='GG = better fat oxidation, longer work capacity',
        snpedia_url='https://www.snpedia.com/index.php/Rs4253778'
    ),

    'rs2016520': AthleticSNP(
        rsid='rs2016520', gene='PPARD', name='PPAR-delta',
        index='durability',
        favorable_alleles=['C'],  # Higher PPARD expression
        unfavorable_alleles=['T'],
        effect_description='Fat oxidation, slow-fiber gene activation',
        evidence='C allele in elite endurance, increased transcription',
        weight=0.1,
        notes='CC = higher PPARD, better fat oxidation, work capacity',
        snpedia_url='https://www.snpedia.com/index.php/Rs2016520'
    ),

    'rs1799983': AthleticSNP(
        rsid='rs1799983', gene='NOS3', name='Glu298Asp',
        index='durability',
        favorable_alleles=['G'],  # Glu (better vascular response)
        unfavorable_alleles=['T'],  # Asp
        effect_description='Vascular responsiveness, training adaptation',
        evidence='Training response rather than baseline performance',
        weight=0.1,
        notes='GG = better vascular training response, work capacity adaptation'
    ),

    'rs8192678': AthleticSNP(
        rsid='rs8192678', gene='PPARGC1A', name='Gly482Ser (PGC-1α)',
        index='durability',
        favorable_alleles=['G'],  # Gly (better mitochondrial)
        unfavorable_alleles=['A'],  # Ser
        effect_description='Mitochondrial biogenesis, oxidative metabolism',
        evidence='VO2max response, work capacity',
        weight=0.15,
        notes='GG = better mitochondrial capacity, work tolerance',
        snpedia_url='https://www.snpedia.com/index.php/Rs8192678'
    ),

    'rs1544410': AthleticSNP(
        rsid='rs1544410', gene='VDR', name='BsmI',
        index='durability',
        favorable_alleles=['A'],  # Better bone density
        unfavorable_alleles=['G'],
        effect_description='Vitamin D receptor, bone density',
        evidence='Additive with FokI for bone/durability',
        weight=0.1,
        notes='Weaker than FokI, but additive for bone density'
    ),

    'rs731236': AthleticSNP(
        rsid='rs731236', gene='VDR', name='TaqI',
        index='durability',
        favorable_alleles=['A'],  # Better bone density
        unfavorable_alleles=['G'],
        effect_description='Vitamin D receptor, bone density',
        evidence='Additive with FokI for bone/durability',
        weight=0.1,
        notes='Weaker than FokI, but additive for bone density'
    ),

    'rs1695': AthleticSNP(
        rsid='rs1695', gene='GSTP1', name='Ile105Val',
        index='durability',
        favorable_alleles=['A'],  # Ile (2-3x more stable)
        unfavorable_alleles=['G'],  # Val
        effect_description='Glutathione S-transferase, oxidative stress clearance',
        evidence='Val105 enzyme 2-3x less stable than Ile105',
        weight=0.1,
        notes='AA (Ile/Ile) = better detox, recovery capacity',
        snpedia_url='https://www.snpedia.com/index.php/Rs1695'
    ),

    'rs1800795': AthleticSNP(
        rsid='rs1800795', gene='IL6', name='-174G>C',
        index='durability',
        favorable_alleles=['C'],  # Lower IL-6
        unfavorable_alleles=['G'],
        effect_description='IL-6 production, inflammation',
        evidence='Inflammation, recovery capacity',
        weight=0.03,  # Reduced: correlated with IL6R/CRP
        notes='CC = lower IL-6, less inflammation'
    ),

    'rs1800629': AthleticSNP(
        rsid='rs1800629', gene='TNF', name='-308G>A',
        index='durability',
        favorable_alleles=['G'],  # Lower TNF
        unfavorable_alleles=['A'],
        effect_description='TNF-alpha production, muscle damage',
        evidence='Inflammation, muscle damage, recovery',
        weight=0.03,  # Reduced: correlated with other inflammation markers
        notes='GG = lower TNF, better recovery from damage'
    ),

    # ============================================================================
    # ADAPTATION INDEX - Learning rate, NOT strength
    # ============================================================================

    'rs6265': AthleticSNP(
        rsid='rs6265', gene='BDNF', name='Val66Met',
        index='adaptation',
        favorable_alleles=['G'],  # Val (better motor learning)
        unfavorable_alleles=['A'],  # Met (impaired learning)
        effect_description='↑ learning: Motor skill acquisition, neuroplasticity',
        evidence='Motor learning studies, NOT endurance/power',
        weight=0.1,
        notes='GG (Val/Val) = better motor learning; affects skill rate, not strength',
        snpedia_url='https://www.snpedia.com/index.php/Rs6265'
    ),

    'rs4680': AthleticSNP(
        rsid='rs4680', gene='COMT', name='Val158Met',
        index='adaptation',
        favorable_alleles=['G'],  # Val (stress resilient, "warrior")
        unfavorable_alleles=['A'],  # Met (stress vulnerable, "worrier")
        effect_description='↑ learning: Stress resilience, dopamine regulation',
        evidence='Warrior/Worrier - affects stress response, weak athletic link',
        weight=0.1,
        notes='GG (Val/Val) = stress resilient; psychological trait, not performance',
        snpedia_url='https://www.snpedia.com/index.php/Rs4680'
    ),

    'rs1800497': AthleticSNP(
        rsid='rs1800497', gene='ANKK1', name='Taq1A (DRD2)',
        index='adaptation',
        favorable_alleles=['C'],  # A2 (normal D2 receptors)
        unfavorable_alleles=['T'],  # A1 (40% reduced receptors)
        effect_description='↑ learning: Motivation, reward processing',
        evidence='Affects motivation for ALL sports equally',
        weight=0.05,
        notes='CC (A2/A2) = better motivation; not endurance/power specific',
        snpedia_url='https://www.snpedia.com/index.php/Rs1800497'
    ),

    # ============================================================================
    # METABOLIC INDEX - Body comp, health-adjacent, explicitly non-performance
    # ============================================================================

    'rs9939609': AthleticSNP(
        rsid='rs9939609', gene='FTO', name='Fat mass and obesity',
        index='metabolic',
        favorable_alleles=['T'],  # Protective, lower BMI
        unfavorable_alleles=['A'],  # Risk allele for obesity
        effect_description='↑ metabolic efficiency: Obesity risk, appetite control',
        evidence='A allele OR=1.34-1.55 for obesity',
        weight=0.2,
        notes='TT = lower obesity risk; health marker, indirect performance',
        snpedia_url='https://www.snpedia.com/index.php/Rs9939609'
    ),

    'rs659366': AthleticSNP(
        rsid='rs659366', gene='UCP2', name='-866G>A',
        index='metabolic',
        favorable_alleles=['A'],  # More efficient
        unfavorable_alleles=['G'],
        effect_description='↑ metabolic efficiency: Mitochondrial uncoupling',
        evidence='Energy expenditure, efficiency',
        weight=0.1,
        notes='AA = more efficient energy use; indirect performance'
    ),

    'rs4994': AthleticSNP(
        rsid='rs4994', gene='ADRB3', name='Trp64Arg',
        index='metabolic',
        favorable_alleles=['T'],  # Trp (better lipolysis)
        unfavorable_alleles=['C'],  # Arg (lower lipolysis)
        effect_description='↑ metabolic efficiency: Lipolysis, fat metabolism',
        evidence='Fat metabolism, body composition',
        weight=0.1,
        notes='TT = better lipolysis; health marker, indirect performance'
    ),

    'rs1800849': AthleticSNP(
        rsid='rs1800849', gene='UCP3', name='UCP3 -55C>T',
        index='metabolic',
        favorable_alleles=['T'],  # Lower BMI tendency
        unfavorable_alleles=['C'],  # Higher BMI
        effect_description='↑ metabolic efficiency: Thermogenesis, energy expenditure',
        evidence='UCP3 promoter, cold adaptation, BMI associations',
        weight=0.1,
        notes='TT = lower BMI in some populations; indirect performance',
        snpedia_url='https://www.snpedia.com/index.php/Rs1800849'
    ),

    'rs470117': AthleticSNP(
        rsid='rs470117', gene='CPT1B', name='E531K',
        index='metabolic',
        favorable_alleles=['A'],  # Normal beta-oxidation
        unfavorable_alleles=['G'],  # Decreased oxidation
        effect_description='↑ metabolic efficiency: Fatty acid transport',
        evidence='G allele decreases mitochondrial beta-oxidation',
        weight=0.1,
        notes='AA = normal beta-oxidation; indirect performance',
        snpedia_url='https://www.snpedia.com/index.php/Rs470117'
    ),

    # Adrenergic - metabolic context only (not performance in non-asthmatics)
    'rs1042713': AthleticSNP(
        rsid='rs1042713', gene='ADRB2', name='Arg16Gly',
        index='metabolic',
        favorable_alleles=['G'],  # Gly
        unfavorable_alleles=['A'],  # Arg
        effect_description='Beta-2 receptor, exercise economy',
        evidence='Inconsistent evidence across studies',
        weight=0.05,
        notes='Weak performance association'
    ),

    'rs1042714': AthleticSNP(
        rsid='rs1042714', gene='ADRB2', name='Gln27Glu',
        index='metabolic',
        favorable_alleles=['G'],  # Glu (resists downregulation)
        unfavorable_alleles=['C'],  # Gln
        effect_description='Beta-2 receptor downregulation resistance',
        evidence='Bronchodilation benefit, weak performance link',
        weight=0.05,
        notes='GG = better bronchodilation; weak athletic association',
        snpedia_url='https://www.snpedia.com/index.php/Rs1042714'
    ),

    'rs1801253': AthleticSNP(
        rsid='rs1801253', gene='ADRB1', name='Arg389Gly',
        index='metabolic',
        favorable_alleles=['C'],  # Arg389 (higher activity)
        unfavorable_alleles=['G'],  # Gly389
        effect_description='Beta-1 receptor, cardiac contractility',
        evidence='Cardiac function, weak performance link outside cardio',
        weight=0.05,
        notes='CC = better cardiac response; adrenergic ≠ performance in non-asthmatics',
        snpedia_url='https://www.snpedia.com/index.php/Rs1801253'
    ),

    'rs1049434': AthleticSNP(
        rsid='rs1049434', gene='SLC16A1', name='MCT1 A1470T',
        index='metabolic',
        favorable_alleles=['T'],  # Better lactate transport
        unfavorable_alleles=['A'],
        effect_description='Lactate transport',
        evidence='Lactate clearance, fatigue resistance',
        weight=0.05,
        notes='TT = better lactate transport'
    ),

    'rs35767': AthleticSNP(
        rsid='rs35767', gene='IGF1', name='IGF-1 promoter',
        index='output',
        favorable_alleles=['T'],  # Raises IGF-1
        unfavorable_alleles=['C'],  # Lower IGF-1
        effect_description='IGF-1 regulation, hypertrophy',
        evidence='T allele raises IGF-1, hypertrophy advantage',
        weight=0.1,
        notes='TT = higher IGF-1, faster hypertrophy',
        snpedia_url='https://www.snpedia.com/index.php/Rs35767'
    ),

    'rs8111989': AthleticSNP(
        rsid='rs8111989', gene='CKM', name='Muscle creatine kinase',
        index='durability',
        favorable_alleles=['A'],  # Less muscle damage
        unfavorable_alleles=['G'],
        effect_description='Creatine kinase, muscle damage/recovery',
        evidence='Marathon performance, preliminary',
        weight=0.05,
        notes='AA = less muscle damage'
    ),

    'rs2070744': AthleticSNP(
        rsid='rs2070744', gene='NOS3', name='-786T>C',
        index='durability',
        favorable_alleles=['T'],  # Higher NO
        unfavorable_alleles=['C'],
        effect_description='Nitric oxide production',
        evidence='Vascular function, weak performance link',
        weight=0.05,
        notes='TT = higher NO production'
    ),

    'rs2010963': AthleticSNP(
        rsid='rs2010963', gene='VEGFA', name='VEGF -634G>C',
        index='durability',
        favorable_alleles=['C'],  # Higher VEGF
        unfavorable_alleles=['G'],
        effect_description='Vascular endothelial growth factor',
        evidence='Capillarization, oxygen delivery',
        weight=0.05,
        notes='CC = better capillarization'
    ),

    # Dropped: HIF1A rs11549465, EPAS1 rs1867785/rs56721780
    # Reason: Ancestry-gated noise, redundant with PPARGC1A, 0.05 weight adds no info

    # Dropped: ADRB1 rs1801253
    # Reason: Moved to metabolic - adrenergic ≠ performance in non-asthmatics
}


class AthleticScorer:
    """Calculate 4 orthogonal athletic indices"""

    def __init__(self, snps_data: snps.SNPs):
        self.snps_data = snps_data
        self.results = {}
        self.missing_snps = []

    @staticmethod
    def complement(allele: str) -> str:
        """Return complement of DNA base"""
        comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return comp.get(allele.upper(), allele)

    def calculate_snp_score(self, snp_info: AthleticSNP) -> Tuple[Optional[float], str, str]:
        """
        Calculate weighted score for a single SNP
        Returns: (score, genotype, interpretation)
        """
        snp_data = self.snps_data.snps[self.snps_data.snps.index == snp_info.rsid]

        if len(snp_data) == 0 or snp_data['genotype'].isna().all():
            return None, 'Missing', 'Data not available'

        genotype = str(snp_data['genotype'].iloc[0])

        # Create complement allele sets for strand flip compatibility
        favorable = set(snp_info.favorable_alleles)
        unfavorable = set(snp_info.unfavorable_alleles)
        favorable_complements = {self.complement(a) for a in favorable}
        unfavorable_complements = {self.complement(a) for a in unfavorable}

        # Count alleles (checking both strands)
        favorable_count = sum(1 for a in genotype if a in favorable or a in favorable_complements)
        unfavorable_count = sum(1 for a in genotype if a in unfavorable or a in unfavorable_complements)

        # Calculate weighted score
        if favorable_count == 2:
            score = snp_info.weight  # Homozygous favorable
        elif unfavorable_count == 2:
            score = -snp_info.weight  # Homozygous unfavorable
        elif favorable_count == 1 and unfavorable_count == 1:
            score = 0  # Heterozygous - neutral
        else:
            score = 0  # Unknown

        # Check strand orientation
        is_flipped = not any(a in favorable | unfavorable for a in genotype)
        if is_flipped:
            strand_info = f"{genotype} [flip: {''.join(self.complement(a) for a in genotype)}]"
        else:
            strand_info = genotype

        # Generate interpretation
        if favorable_count == 2:
            interpretation = f"Favorable ({strand_info})"
        elif unfavorable_count == 2:
            interpretation = f"Unfavorable ({strand_info})"
        elif favorable_count == 1:
            interpretation = f"Neutral/het ({strand_info})"
        else:
            interpretation = f"Neutral ({strand_info})"

        return score, genotype, interpretation

    def calculate_indices(self) -> Dict:
        """Calculate the 4 orthogonal indices"""

        index_scores = {
            'output': 0,
            'durability': 0,
            'adaptation': 0,
            'metabolic': 0
        }

        snp_scores = []

        for rsid, snp_info in ATHLETIC_SNPS.items():
            score, genotype, interpretation = self.calculate_snp_score(snp_info)

            if score is not None:
                index_scores[snp_info.index] += score

                snp_scores.append({
                    'rsid': rsid,
                    'gene': snp_info.gene,
                    'name': snp_info.name,
                    'index': snp_info.index,
                    'score': score,
                    'weight': snp_info.weight,
                    'genotype': genotype,
                    'interpretation': interpretation,
                    'effect': snp_info.effect_description,
                    'notes': snp_info.notes
                })
            else:
                self.missing_snps.append(rsid)

        return {
            'index_scores': index_scores,
            'snp_scores': snp_scores,
            'missing_snps': self.missing_snps,
            'snps_analyzed': len(snp_scores),
            'snps_total': len(ATHLETIC_SNPS)
        }

    def print_report(self, results: Dict):
        """Print formatted 4-index report"""

        # Color codes
        HEADER = '\033[95m'
        BLUE = '\033[94m'
        CYAN = '\033[96m'
        GREEN = '\033[92m'
        YELLOW = '\033[93m'
        RED = '\033[91m'
        BOLD = '\033[1m'
        END = '\033[0m'

        print(f"\n{BOLD}{CYAN}{'='*80}{END}")
        print(f"{BOLD}{CYAN} TRAINING BIAS & DURABILITY PROFILE (LIMITED SNP PANEL){END}")
        print(f"{BOLD}{CYAN}{'='*80}{END}\n")

        print(f"{BOLD}SNPs Analyzed:{END} {results['snps_analyzed']}/{results['snps_total']}")

        # 4 Orthogonal Indices
        print(f"\n{BOLD}{YELLOW}═══ ORTHOGONAL INDICES ═══{END}\n")

        indices = results['index_scores']

        # OUTPUT INDEX
        output_score = indices['output']
        output_color = GREEN if output_score > 0.5 else RED if output_score < -0.5 else YELLOW
        print(f"{BOLD}OUTPUT INDEX:{END} {output_color}{output_score:+.2f}{END}")
        print(f"  Force/speed/power expression (what wins races/lifts)")
        if output_score > 0.8:
            print(f"  → {GREEN}Strong power potential - explosive advantage{END}")
        elif output_score > 0.3:
            print(f"  → {GREEN}Moderate power advantage{END}")
        elif output_score < -0.8:
            print(f"  → {RED}Power disadvantage - train technique & leverage{END}")
        elif output_score < -0.3:
            print(f"  → {YELLOW}Mild power disadvantage{END}")
        else:
            print(f"  → {YELLOW}Neutral - trainable either way{END}")

        # DURABILITY INDEX
        durability_score = indices['durability']
        durability_color = GREEN if durability_score > 0.5 else RED if durability_score < -0.5 else YELLOW
        print(f"\n{BOLD}DURABILITY INDEX:{END} {durability_color}{durability_score:+.2f}{END}")
        print(f"  Work capacity, injury resistance (training tolerance)")
        if durability_score > 1.0:
            print(f"  → {GREEN}Excellent work capacity - high training volume tolerated{END}")
        elif durability_score > 0.3:
            print(f"  → {GREEN}Good durability - above-average training tolerance{END}")
        elif durability_score < -1.0:
            print(f"  → {RED}Fragility risk - prioritize recovery, lower volume{END}")
        elif durability_score < -0.3:
            print(f"  → {YELLOW}Below-average durability - careful progression{END}")
        else:
            print(f"  → {YELLOW}Average durability{END}")

        # ADAPTATION INDEX
        adaptation_score = indices['adaptation']
        adaptation_color = GREEN if adaptation_score > 0.1 else RED if adaptation_score < -0.1 else YELLOW
        print(f"\n{BOLD}ADAPTATION INDEX:{END} {adaptation_color}{adaptation_score:+.2f}{END}")
        print(f"  Skill acquisition, neuroplasticity (learning rate, NOT strength)")
        if adaptation_score > 0.15:
            print(f"  → {GREEN}Fast learner - technical skills improve quickly{END}")
        elif adaptation_score < -0.15:
            print(f"  → {RED}Slower learning - needs more reps for skill mastery{END}")
        else:
            print(f"  → {YELLOW}Average learning rate{END}")

        # METABOLIC INDEX
        metabolic_score = indices['metabolic']
        metabolic_color = GREEN if metabolic_score > 0.2 else RED if metabolic_score < -0.2 else YELLOW
        print(f"\n{BOLD}METABOLIC INDEX:{END} {metabolic_color}{metabolic_score:+.2f}{END}")
        print(f"  Body composition, energy efficiency (health-adjacent, indirect)")
        if metabolic_score > 0.3:
            print(f"  → {GREEN}Favorable body comp genetics - easier to stay lean{END}")
        elif metabolic_score < -0.3:
            print(f"  → {RED}Obesity-prone - diet/activity critical{END}")
        else:
            print(f"  → {YELLOW}Average metabolic genetics{END}")

        # Detailed SNP results by index
        print(f"\n{BOLD}{YELLOW}═══ DETAILED SNP RESULTS ═══{END}\n")

        # Special handling for durability - group by mechanism
        durability_groups = {
            'Connective Tissue': ['COL5A1', 'COL1A1', 'COL1A2', 'COL12A1', 'COL3A1'],
            'Oxidative Stress & Recovery': ['SOD2', 'NFE2L2', 'GSTP1'],
            'Inflammation': ['IL6R', 'CRP', 'IL6', 'TNF'],
            'Metabolic Capacity': ['PPARA', 'PPARD', 'PPARGC1A', 'NOS3'],
            'Other': []  # catch-all
        }

        for index_name in ['output', 'durability', 'adaptation', 'metabolic']:
            index_snps = [s for s in results['snp_scores'] if s['index'] == index_name]
            if index_snps:
                index_label = index_name.upper() + " INDEX"
                print(f"{BOLD}{CYAN}{index_label}:{END}")

                if index_name == 'durability':
                    # Group durability SNPs by mechanism
                    for group_name, gene_list in durability_groups.items():
                        group_snps = [s for s in index_snps if s['gene'] in gene_list]
                        if group_snps:
                            print(f"  {BOLD}{group_name}:{END}")
                            for snp in group_snps:
                                score_symbol = '→' if snp['score'] == 0 else '▲' if snp['score'] > 0 else '▼'
                                score_color = GREEN if snp['score'] > 0 else RED if snp['score'] < 0 else YELLOW
                                print(f"    {snp['gene']} {snp['name']} ({snp['rsid']}) [wt:{snp['weight']:.2f}] {score_color}{score_symbol} {snp['score']:+.2f}{END}")
                                # Show SNPedia link for high-weight SNPs (≥0.15)
                                snp_obj = ATHLETIC_SNPS.get(snp['rsid'])
                                if snp['weight'] >= 0.15 and snp_obj and snp_obj.snpedia_url:
                                    print(f"      SNPedia: {snp_obj.snpedia_url}")
                                if snp['notes']:
                                    print(f"      {snp['notes']}")

                    # Catch remaining durability SNPs not in groups
                    remaining = [s for s in index_snps if s['gene'] not in sum(durability_groups.values(), [])]
                    if remaining:
                        print(f"  {BOLD}Other:{END}")
                        for snp in remaining:
                            score_symbol = '→' if snp['score'] == 0 else '▲' if snp['score'] > 0 else '▼'
                            score_color = GREEN if snp['score'] > 0 else RED if snp['score'] < 0 else YELLOW
                            print(f"    {snp['gene']} {snp['name']} ({snp['rsid']}) [wt:{snp['weight']:.2f}] {score_color}{score_symbol} {snp['score']:+.2f}{END}")
                            if snp['notes']:
                                print(f"      {snp['notes']}")
                else:
                    # Normal display for other indices
                    for snp in index_snps:
                        score_symbol = '→' if snp['score'] == 0 else '▲' if snp['score'] > 0 else '▼'
                        score_color = GREEN if snp['score'] > 0 else RED if snp['score'] < 0 else YELLOW

                        print(f"  {BOLD}{snp['gene']}{END} {snp['name']} ({snp['rsid']}) [weight: {snp['weight']:.2f}]")
                        print(f"    {score_color}{score_symbol} {snp['interpretation']}{END} (score: {snp['score']:+.2f})")
                        # Show SNPedia link for high-weight SNPs (≥0.15)
                        snp_obj = ATHLETIC_SNPS.get(snp['rsid'])
                        if snp['weight'] >= 0.15 and snp_obj and snp_obj.snpedia_url:
                            print(f"    SNPedia: {snp_obj.snpedia_url}")
                        if snp['notes']:
                            print(f"    {snp['notes']}")
                print()

        # Missing SNPs
        if results['missing_snps']:
            print(f"{BOLD}{YELLOW}Missing SNPs:{END}")
            print(f"  {', '.join(results['missing_snps'])[:200]}...")
            print()

        print(f"{BOLD}{RED}⚠ CRITICAL LIMITATIONS:{END}")
        print(f"{YELLOW}• This analyzes {results['snps_analyzed']} SNPs from a limited panel{END}")
        print(f"{YELLOW}• Missing hundreds of GWAS-identified variants{END}")
        print(f"{YELLOW}• NO PHENOTYPE ANCHORING - genetics without training data is astrology{END}")
        print(f"{YELLOW}• Each index measures ONE causal mechanism - don't average them{END}")
        print(f"{YELLOW}• Genetics explains ~20-30% of athletic variance{END}")
        print(f"\n{BOLD}Add phenotype data (lifts, times, injury history) for Bayesian updating.{END}")
        print(f"{BOLD}Training, technique, and psychology >>> genetics.{END}\n")


def analyze_athletic_performance(snps_data: snps.SNPs) -> Dict:
    """Main entry point for athletic analysis"""
    scorer = AthleticScorer(snps_data)
    results = scorer.calculate_indices()
    scorer.print_report(results)
    return results


if __name__ == "__main__":
    import sys

    filepath = "AncestryDNA.txt"
    if len(sys.argv) > 1:
        filepath = sys.argv[1]

    try:
        s = snps.SNPs(filepath)
        results = analyze_athletic_performance(s)
    except FileNotFoundError:
        print(f"Error: Could not find file '{filepath}'")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
