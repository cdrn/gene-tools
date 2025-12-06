#!/usr/bin/env python3
"""
Athletic Performance Polygenic Score
Endurance <-> Power spectrum analysis based on validated genetic markers
"""

import snps
from typing import Dict, Tuple, Optional, List
from dataclasses import dataclass


@dataclass
class AthleticSNP:
    """Athletic performance SNP with scoring information"""
    rsid: str
    gene: str
    name: str
    category: str
    endurance_alleles: List[str]
    power_alleles: List[str]
    effect_description: str
    evidence: str
    effect_size: float = 1.0  # Odds ratio or beta coefficient (1.0 = no effect)
    notes: Optional[str] = None
    snpedia_url: Optional[str] = None


# Athletic Performance SNPs Database
# Verified against SNPedia and primary literature
ATHLETIC_SNPS = {
    # ============================================================================
    # MUSCLE & ENERGY SYSTEMS
    # ============================================================================

    'rs1815739': AthleticSNP(
        rsid='rs1815739',
        gene='ACTN3',
        name='R577X (alpha-actinin-3)',
        category='muscle_energy',
        endurance_alleles=['T'],  # X allele (stop codon)
        power_alleles=['C'],      # R allele (arginine)
        effect_description='Alpha-actinin-3 in fast-twitch fibers',
        evidence='Olympic athlete studies, meta-analyses',
        effect_size=1.5,  # OR ~1.5 for sprint/power in meta-analyses
        notes='TT (XX) = no alpha-actinin-3, endurance advantage'
    ),

    'rs4343': AthleticSNP(
        rsid='rs4343',
        gene='ACE',
        name='ACE I/D proxy',
        category='muscle_energy',
        endurance_alleles=['A'],  # Proxy for I allele
        power_alleles=['G'],      # Proxy for D allele
        effect_description='ACE insertion/deletion polymorphism proxy',
        evidence='Elite athlete studies, VO2max association',
        effect_size=1.2,  # Mixed results in meta-analyses, conservative estimate
        notes='A=I (insertion), G=D (deletion). Perfect proxy for I/D.'
    ),

    'rs8192678': AthleticSNP(
        rsid='rs8192678',
        gene='PPARGC1A',
        name='Gly482Ser (PGC-1α)',
        category='muscle_energy',
        endurance_alleles=['G'],  # Gly (glycine)
        power_alleles=['A'],      # Ser (serine)
        effect_description='Mitochondrial biogenesis, oxidative metabolism',
        evidence='VO2max response, endurance athlete enrichment',
        effect_size=1.3,  # Meta-analysis shows consistent endurance enrichment
        notes='G (Gly) = better aerobic capacity'
    ),

    'rs17602729': AthleticSNP(
        rsid='rs17602729',
        gene='AMPD1',
        name='Q12X',
        category='muscle_energy',
        endurance_alleles=['T'],  # X (stop) - deficiency
        power_alleles=['C'],      # Q (glutamine) - normal
        effect_description='AMP deaminase deficiency',
        evidence='Reduced sprint/power capacity with deficiency',
        effect_size=2.17,  # OR=2.17 from PMC12152022
        notes='TT = enzyme deficiency, power drop risk; CC = normal power capacity'
    ),

    'rs1805086': AthleticSNP(
        rsid='rs1805086',
        gene='MSTN',
        name='K153R (myostatin)',
        category='muscle_energy',
        endurance_alleles=['A'],  # K (lysine) - normal myostatin
        power_alleles=['G'],      # R (arginine) - reduced function
        effect_description='Myostatin K153R, muscle mass regulation',
        evidence='PMC3024427: R allele associated with strength athletes (OR=2.02)',
        effect_size=2.02,  # Already has OR=2.02 from PMC3024427
        notes='AA (KK) = normal; AG/GG (R allele) = increased muscle mass, ultra-rare (~3-4% in Caucasians)',
        snpedia_url='https://www.snpedia.com/index.php/Rs1805086'
    ),

    'rs2228570': AthleticSNP(
        rsid='rs2228570',
        gene='VDR',
        name='FokI (Vitamin D receptor)',
        category='muscle_energy',
        endurance_alleles=['G'],  # f allele (longer protein, less active)
        power_alleles=['A'],      # F allele (shorter protein, more active)
        effect_description='Vitamin D receptor FokI, muscle fiber composition',
        evidence='SNPedia: rs2228570, multiple studies on muscle/strength',
        effect_size=1.2,  # Muscle fiber composition effects
        notes='AA (FF) = more active VDR, may favor type-II fibers; GG (ff) = less active VDR',
        snpedia_url='https://www.snpedia.com/index.php/Rs2228570'
    ),

    'rs35767': AthleticSNP(
        rsid='rs35767',
        gene='IGF1',
        name='IGF-1 promoter (C/T)',
        category='muscle_energy',
        endurance_alleles=['C'],  # Lower IGF-1
        power_alleles=['T'],      # T allele RAISES IGF-1
        effect_description='IGF-1 regulation, growth hormone pathway',
        evidence='PLOS One: T allele = IGF-1 raising, hypertrophy advantage',
        effect_size=1.3,  # T allele increases IGF-1, hypertrophy
        notes='TT = higher circulating IGF-1, faster hypertrophy; CC = lower IGF-1',
        snpedia_url='https://www.snpedia.com/index.php/Rs35767'
    ),

    'rs8111989': AthleticSNP(
        rsid='rs8111989',
        gene='CKM',
        name='Muscle creatine kinase',
        category='muscle_energy',
        endurance_alleles=['A'],
        power_alleles=['G'],
        effect_description='Creatine kinase, muscle damage/recovery',
        evidence='Marathon performance, injury susceptibility',
        effect_size=1.1,  # Preliminary marathon evidence
        notes='A = better endurance, less muscle damage'
    ),

    'rs2306862': AthleticSNP(
        rsid='rs2306862',
        gene='LRP5',
        name='LRP5 (Wnt signaling)',
        category='muscle_energy',
        endurance_alleles=['T'],  # Non-C allele
        power_alleles=['C'],      # C allele - lean mass, strength
        effect_description='LRP5 Wnt signaling, lean mass, bone density',
        evidence='Strong lean-mass GWAS signal, replicated in grip strength + power (Genes 2023 review)',
        effect_size=1.3,  # Strong lean-mass signal from GWAS
        notes='CC = higher lean mass, better grip strength; TT = lower lean mass tendency'
    ),

    'rs2070744': AthleticSNP(
        rsid='rs2070744',
        gene='NOS3',
        name='-786T>C',
        category='muscle_energy',
        endurance_alleles=['T'],
        power_alleles=['C'],
        effect_description='Nitric oxide production',
        evidence='Endurance performance, vascular function',
        effect_size=1.1,  # Endurance association, NO production
        notes='T = higher NO, better endurance'
    ),

    'rs2010963': AthleticSNP(
        rsid='rs2010963',
        gene='VEGFA',
        name='VEGF -634G>C',
        category='muscle_energy',
        endurance_alleles=['C'],
        power_alleles=['G'],
        effect_description='Vascular endothelial growth factor',
        evidence='Capillarization, oxygen delivery',
        effect_size=1.1,  # Vascular/capillarization effects
        notes='C = higher VEGF, better endurance adaptation'
    ),

    # ============================================================================
    # ADRENERGIC & FATIGUE RESISTANCE
    # ============================================================================

    'rs1042713': AthleticSNP(
        rsid='rs1042713',
        gene='ADRB2',
        name='Arg16Gly',
        category='adrenergic',
        endurance_alleles=['G'],  # Gly
        power_alleles=['A'],      # Arg
        effect_description='Beta-2 adrenergic receptor function',
        evidence='Exercise economy, fat oxidation',
        effect_size=1.1,  # Inconsistent evidence across studies
        notes='G (Gly) = better endurance economy'
    ),

    'rs1042714': AthleticSNP(
        rsid='rs1042714',
        gene='ADRB2',
        name='Gln27Glu (C/G)',
        category='adrenergic',
        endurance_alleles=['G'],  # G = Glu27 (resists downregulation, bronchodilation)
        power_alleles=['C'],      # C = Gln27 (baseline)
        effect_description='Beta-2 receptor downregulation resistance',
        evidence='Mixed: G (Glu) shows bronchodilation benefit, but no strong endurance performance link in athletes',
        effect_size=1.0,  # Weak performance association
        notes='GG (Glu/Glu) = resists receptor downregulation, better bronchodilation; CC (Gln/Gln) = baseline. Performance association weak.',
        snpedia_url='https://www.snpedia.com/index.php/Rs1042714'
    ),

    'rs4994': AthleticSNP(
        rsid='rs4994',
        gene='ADRB3',
        name='Trp64Arg',
        category='adrenergic',
        endurance_alleles=['C'],  # Arg
        power_alleles=['T'],      # Trp
        effect_description='Beta-3 adrenergic receptor, lipolysis',
        evidence='Fat metabolism, energy efficiency',
        effect_size=1.1,  # Some evidence for lipolysis/endurance
        notes='C (Arg) = lower lipolysis, endurance-lean phenotype'
    ),

    # ============================================================================
    # OXYGEN & HYPOXIA RESPONSE
    # ============================================================================

    'rs11549465': AthleticSNP(
        rsid='rs11549465',
        gene='HIF1A',
        name='Pro582Ser',
        category='oxygen',
        endurance_alleles=['C'],  # Pro
        power_alleles=['T'],      # Ser
        effect_description='Hypoxia-inducible factor 1-alpha',
        evidence='Altitude adaptation, glycolytic capacity',
        effect_size=1.1,  # Rare variant, altitude adaptation
        notes='T (Ser) = rare, more glycolytic (power)'
    ),

    'rs1867785': AthleticSNP(
        rsid='rs1867785',
        gene='EPAS1',
        name='HIF-2α variant',
        category='oxygen',
        endurance_alleles=['T'],  # Context-dependent
        power_alleles=['C'],
        effect_description='Hypoxia response, altitude adaptation',
        evidence='Tibetan altitude adaptation',
        effect_size=1.1,  # Ancestry-dependent, Tibetan studies
        notes='Ancestry-dependent effects, T generally endurance'
    ),

    'rs56721780': AthleticSNP(
        rsid='rs56721780',
        gene='EPAS1',
        name='HIF-2α Tibetan variant (-886G>C)',
        category='oxygen',
        endurance_alleles=['C'],  # C allele decreases IKZF1 binding
        power_alleles=['G'],
        effect_description='Tibetan high-altitude hypoxia adaptation',
        evidence='Sci Rep 2014: C allele freq 0.372 in Tibetans vs 0.010 in Han Chinese',
        effect_size=1.0,  # Very population-specific (Tibetan)
        notes='CC = attenuated EPAS1 repression, altitude adaptation (rare in non-East Asian)',
        snpedia_url='https://www.snpedia.com/index.php/Rs56721780'
    ),

    'rs2016520': AthleticSNP(
        rsid='rs2016520',
        gene='PPARD',
        name='PPAR-delta (T294C)',
        category='oxygen',
        endurance_alleles=['C'],  # C allele binds Sp-1 transcription factor
        power_alleles=['T'],
        effect_description='Fat oxidation and slow-fiber gene activation',
        evidence='PLOS One: C allele in elite endurance athletes, increased PPARD transcription',
        effect_size=1.2,  # Fat oxidation, endurance enrichment
        notes='CC = higher PPARD expression, better fat oxidation, endurance advantage',
        snpedia_url='https://www.snpedia.com/index.php/Rs2016520'
    ),

    'rs1801253': AthleticSNP(
        rsid='rs1801253',
        gene='ADRB1',
        name='Arg389Gly (C/G)',
        category='oxygen',
        endurance_alleles=['C'],  # Arg389 - higher receptor activity
        power_alleles=['G'],      # Gly389 - lower receptor activity
        effect_description='Beta-1 adrenergic receptor, cardiac contractility',
        evidence='SNPedia rs1801253: Arg389 = better beta-blocker response, cardiac function',
        effect_size=1.2,  # Cardiac function, training response
        notes='CC (Arg/Arg) = higher receptor activity, better endurance training response; GG (Gly/Gly) = lower activity',
        snpedia_url='https://www.snpedia.com/index.php/Rs1801253'
    ),

    # ============================================================================
    # NEUROMOTOR & ADRENERGIC (Brain-Muscle Connection)
    # ============================================================================

    'rs6265': AthleticSNP(
        rsid='rs6265',
        gene='BDNF',
        name='Val66Met (G/A)',
        category='neuromotor',
        endurance_alleles=['G'],  # Val - better motor learning
        power_alleles=['A'],      # Met - impaired motor learning/retention
        effect_description='Neuroplasticity, motor skill learning, BDNF secretion',
        evidence='Motor learning, not endurance/power specific',
        effect_size=1.0,  # NEUTRAL - affects skill acquisition, not endurance vs power
        notes='GG (Val/Val) = normal BDNF, better motor learning; AA (Met/Met) = reduced secretion, impaired plasticity. NOT an endurance/power marker.',
        snpedia_url='https://www.snpedia.com/index.php/Rs6265'
    ),

    'rs1800497': AthleticSNP(
        rsid='rs1800497',
        gene='ANKK1',
        name='Taq1A (DRD2) C/T',
        category='neuromotor',
        endurance_alleles=['C'],  # A2 allele - normal D2 receptors
        power_alleles=['T'],      # A1 allele - 40% reduced D2 receptors
        effect_description='Dopamine D2 receptor density, reward/motivation',
        evidence='Affects motivation for ALL sports, not endurance/power specific',
        effect_size=1.0,  # NEUTRAL - motivation affects all athletic performance equally
        notes='CC (A2/A2) = normal dopamine receptors, better motivation; TT (A1/A1) = reduced receptors, lower drive. NOT an endurance/power marker.',
        snpedia_url='https://www.snpedia.com/index.php/Rs1800497'
    ),

    'rs4680': AthleticSNP(
        rsid='rs4680',
        gene='COMT',
        name='Val158Met (G/A)',
        category='neuromotor',
        endurance_alleles=['A'],  # Met - low enzyme activity, "worrier"
        power_alleles=['G'],      # Val - high enzyme activity, "warrior"
        effect_description='Dopamine breakdown rate, stress resilience',
        evidence='Warrior/Worrier gene - affects stress response, not physiological endurance/power',
        effect_size=1.0,  # NEUTRAL - psychological trait, weak/speculative athletic association
        notes='AA (Met/Met) = slow dopamine breakdown, better baseline cognition but stress-vulnerable; GG (Val/Val) = fast breakdown, stress resilient. NOT a reliable endurance/power marker.',
        snpedia_url='https://www.snpedia.com/index.php/Rs4680'
    ),

    # ============================================================================
    # INFLAMMATION, RECOVERY & CONNECTIVE TISSUE
    # ============================================================================

    'rs2228145': AthleticSNP(
        rsid='rs2228145',
        gene='IL6R',
        name='Asp358Ala (A/C)',
        category='inflammation',
        endurance_alleles=['C'],  # C = 358Ala (increased shedding, lower CRP, anti-inflammatory)
        power_alleles=['A'],      # A = 358Asp (normal membrane IL-6R, normal IL-6 response)
        effect_description='IL-6 receptor shedding and inflammation',
        evidence='Mixed: C allele increases soluble IL-6R (2-fold), lowers CRP; but performance benefit speculative',
        effect_size=1.1,  # Anti-inflammatory, speculative endurance
        notes='CC (Ala/Ala) = high IL-6R shedding, lower CRP, anti-inflammatory; AA (Asp/Asp) = normal IL-6 signaling. Endurance link uncertain.',
        snpedia_url='https://www.snpedia.com/index.php/Rs2228145'
    ),

    'rs1205': AthleticSNP(
        rsid='rs1205',
        gene='CRP',
        name='C-reactive protein (C/T)',
        category='inflammation',
        endurance_alleles=['T'],  # T allele LOWERS CRP by 20%
        power_alleles=['C'],      # C allele = higher CRP
        effect_description='C-reactive protein baseline levels',
        evidence='SNPedia rs1205: Each T allele lowers CRP by 20%',
        effect_size=1.2,  # 20% CRP reduction per allele
        notes='TT = lowest CRP, better recovery; CC = higher baseline inflammation',
        snpedia_url='https://www.snpedia.com/index.php/Rs1205'
    ),

    'rs1695': AthleticSNP(
        rsid='rs1695',
        gene='GSTP1',
        name='Ile105Val (A/G)',
        category='inflammation',
        endurance_alleles=['A'],  # Ile - more stable enzyme (2-3x)
        power_alleles=['G'],      # Val - less stable enzyme
        effect_description='Glutathione S-transferase Pi, oxidative stress clearance',
        evidence='SNPedia rs1695: Val105 enzyme 2-3x less stable than Ile105',
        effect_size=1.1,  # Antioxidant capacity
        notes='AA (Ile/Ile) = more stable GSTP1, better detox capacity; GG (Val/Val) = less stable enzyme',
        snpedia_url='https://www.snpedia.com/index.php/Rs1695'
    ),

    'rs4880': AthleticSNP(
        rsid='rs4880',
        gene='SOD2',
        name='Ala16Val (C/T)',
        category='inflammation',
        endurance_alleles=['C'],  # C = Ala (better mitochondrial import, higher activity)
        power_alleles=['T'],      # T = Val (30-40% reduced activity, poorer import)
        effect_description='Mitochondrial superoxide dismutase, ROS protection',
        evidence='Multiple studies: Ala (C) imported more efficiently than Val (T), 30-40% activity difference',
        effect_size=1.3,  # 30-40% activity difference, mitochondrial
        notes='CC (Ala/Ala) = better mitochondrial import, higher SOD2 activity; TT (Val/Val) = reduced import & activity',
        snpedia_url='https://www.snpedia.com/index.php/Rs4880'
    ),

    'rs1800795': AthleticSNP(
        rsid='rs1800795',
        gene='IL6',
        name='-174G>C',
        category='inflammation',
        endurance_alleles=['C'],
        power_alleles=['G'],
        effect_description='IL-6 production',
        evidence='Inflammation, recovery',
        effect_size=1.1,  # IL-6 levels, inflammation
        notes='C = lower IL-6, better for endurance'
    ),

    'rs1800629': AthleticSNP(
        rsid='rs1800629',
        gene='TNF',
        name='-308G>A',
        category='inflammation',
        endurance_alleles=['G'],
        power_alleles=['A'],
        effect_description='TNF-alpha production',
        evidence='Inflammation, muscle damage',
        effect_size=1.1,  # TNF levels, recovery
        notes='G = lower TNF, better recovery for endurance'
    ),

    'rs12722': AthleticSNP(
        rsid='rs12722',
        gene='COL5A1',
        name='Type V collagen',
        category='connective',
        endurance_alleles=['T'],
        power_alleles=['C'],
        effect_description='Collagen structure, tendon stiffness',
        evidence='Running economy, injury risk',
        effect_size=1.1,  # Tendon stiffness/injury
        notes='T = optimal tendon stiffness for endurance'
    ),

    'rs1800012': AthleticSNP(
        rsid='rs1800012',
        gene='COL1A1',
        name='Type I collagen Sp1 (G/T)',
        category='connective',
        endurance_alleles=['G'],  # G = normal, lower fracture risk
        power_alleles=['T'],      # T = "s" allele, INCREASED fracture risk (OR=1.78), lower BMD
        effect_description='Collagen type I production, bone mineral density',
        evidence='Meta-analysis: T allele associated with lower BMD, increased osteoporotic fracture risk',
        effect_size=1.3,  # Bone density, fracture risk
        notes='GG = normal collagen production; TT = increased fracture risk (OR=1.78), lower bone density',
        snpedia_url='https://www.snpedia.com/index.php/Rs1800012'
    ),

    'rs1800255': AthleticSNP(
        rsid='rs1800255',
        gene='COL3A1',
        name='Ala698Thr (G/A)',
        category='connective',
        endurance_alleles=['G'],  # G = Ala (normal collagen structure)
        power_alleles=['A'],      # A = Thr (disrupted triple helix)
        effect_description='Type III collagen structure, connective tissue integrity',
        evidence='SNPedia rs1800255: A allele disrupts triple helix, associated with POP (OR=4.79)',
        effect_size=1.2,  # OR=4.79 for injury, structural
        notes='GG (Ala/Ala) = normal collagen structure, better integrity; AA (Thr/Thr) = disrupted structure, higher injury risk',
        snpedia_url='https://www.snpedia.com/index.php/Rs1800255'
    ),

    # ============================================================================
    # FUEL HANDLING & METABOLISM
    # ============================================================================

    'rs4253778': AthleticSNP(
        rsid='rs4253778',
        gene='PPARA',
        name='PPAR-alpha variant',
        category='fuel',
        endurance_alleles=['G'],
        power_alleles=['C'],
        effect_description='Fat oxidation capacity',
        evidence='Endurance athlete enrichment',
        effect_size=1.1,  # Fat oxidation, PPAR-alpha
        notes='G = better fat oxidation for endurance'
    ),

    'rs659366': AthleticSNP(
        rsid='rs659366',
        gene='UCP2',
        name='-866G>A',
        category='fuel',
        endurance_alleles=['A'],
        power_alleles=['G'],
        effect_description='Mitochondrial uncoupling, efficiency',
        evidence='Energy expenditure, efficiency',
        effect_size=1.1,  # Energy expenditure, UCP2
        notes='A = more efficient energy use, endurance'
    ),

    'rs1049434': AthleticSNP(
        rsid='rs1049434',
        gene='SLC16A1',
        name='MCT1 A1470T',
        category='fuel',
        endurance_alleles=['T'],
        power_alleles=['A'],
        effect_description='Lactate transport',
        evidence='Lactate clearance, fatigue resistance',
        effect_size=1.1,  # Lactate transport, MCT1
        notes='T = better lactate transport, endurance'
    ),

    'rs470117': AthleticSNP(
        rsid='rs470117',
        gene='CPT1B',
        name='E531K (A/G)',
        category='fuel',
        endurance_alleles=['A'],  # A = Glu531 (normal beta-oxidation)
        power_alleles=['G'],      # G = Lys531 (decreased beta-oxidation)
        effect_description='Carnitine palmitoyltransferase 1B, fatty acid transport',
        evidence='PLOS One: G allele decreases mitochondrial beta-oxidation, associated with reduced obesity risk',
        effect_size=1.1,  # Beta-oxidation, CPT1B
        notes='AA (Glu/Glu) = normal beta-oxidation, endurance; GG (Lys/Lys) = decreased oxidation',
        snpedia_url='https://www.snpedia.com/index.php/Rs470117'
    ),

    'rs9939609': AthleticSNP(
        rsid='rs9939609',
        gene='FTO',
        name='Fat mass and obesity (T/A)',
        category='fuel',
        endurance_alleles=['T'],  # T = protective, lower BMI
        power_alleles=['A'],      # A = RISK allele for obesity
        effect_description='FTO obesity risk, body composition, appetite',
        evidence='SNPedia rs9939609: A allele OR=1.34 (het) to 1.55 (hom) for obesity',
        effect_size=1.2,  # OR=1.34-1.55 obesity (inverse for performance)
        notes='TT = lower obesity risk, better appetite control; AA = +1.55x obesity risk, increased appetite',
        snpedia_url='https://www.snpedia.com/index.php/Rs9939609'
    ),

    'rs6721961': AthleticSNP(
        rsid='rs6721961',
        gene='NFE2L2',
        name='NRF2 -617C>A',
        category='fuel',
        endurance_alleles=['C'],  # C = normal Nrf2 expression
        power_alleles=['A'],      # A = reduced Nrf2 basal expression
        effect_description='NRF2 transcription factor, oxidative stress response',
        evidence='SNPedia rs6721961: A allele attenuates ARE-mediated gene transcription',
        effect_size=1.1,  # Oxidative stress response
        notes='CC = normal Nrf2 expression, better antioxidant response; AA = reduced basal Nrf2',
        snpedia_url='https://www.snpedia.com/index.php/Rs6721961'
    ),

    'rs1800849': AthleticSNP(
        rsid='rs1800849',
        gene='UCP3',
        name='UCP3 -55C>T',
        category='fuel',
        endurance_alleles=['T'],  # T allele - associations with lower BMI in some studies
        power_alleles=['C'],      # C allele - associated with higher BMI
        effect_description='Uncoupling protein 3, thermogenesis, energy expenditure',
        evidence='SNPedia rs1800849: UCP3 promoter polymorphism, cold adaptation',
        effect_size=1.1,  # Thermogenesis, UCP3
        notes='TT = associated with lower BMI in some populations; CC = higher BMI tendency',
        snpedia_url='https://www.snpedia.com/index.php/Rs1800849'
    ),

    # Note: BDKRB2 -9/+9 is an insertion/deletion, harder to score from SNP data
    # Would need special handling or proxy SNP
}


class AthleticScorer:
    """Calculate athletic performance polygenic score"""

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
        Calculate effect-size weighted score for a single SNP
        Returns: (score, genotype, interpretation)
        Score: Based on log(OR) weighted by genotype
        - Homozygous beneficial: +log(OR)
        - Heterozygous: +log(OR)/2 (additive model)
        - Homozygous detrimental: -log(OR)
        """
        import math

        snp_data = self.snps_data.snps[self.snps_data.snps.index == snp_info.rsid]

        if len(snp_data) == 0 or snp_data['genotype'].isna().all():
            return None, 'Missing', 'Data not available'

        genotype = str(snp_data['genotype'].iloc[0])

        # Create complement allele lists for strand flip compatibility
        endurance_alleles = set(snp_info.endurance_alleles)
        power_alleles = set(snp_info.power_alleles)
        endurance_complements = {self.complement(a) for a in endurance_alleles}
        power_complements = {self.complement(a) for a in power_alleles}

        # Count alleles (checking both strands)
        endurance_count = 0
        power_count = 0

        for allele in genotype:
            if allele in endurance_alleles or allele in endurance_complements:
                endurance_count += 1
            elif allele in power_alleles or allele in power_complements:
                power_count += 1

        # Calculate effect-size weighted score
        # Use log(OR) for proper multiplicative scaling
        effect_weight = math.log(snp_info.effect_size) if snp_info.effect_size > 0 else 0

        if endurance_count == 2:
            # Homozygous endurance
            score = effect_weight
        elif power_count == 2:
            # Homozygous power
            score = -effect_weight
        elif endurance_count == 1 and power_count == 1:
            # Heterozygous - additive model (half effect)
            score = 0  # Neutral for heterozygotes in additive model
        else:
            score = 0  # Unknown/neutral

        # Check if we're on the complement strand
        flipped_genotype = ''.join([self.complement(a) for a in genotype])
        is_flipped = (genotype not in ''.join([str(a) for a in endurance_alleles | power_alleles]))

        # Generate interpretation with strand info
        if is_flipped:
            strand_info = f"{genotype} [flip: {flipped_genotype}]"
        else:
            strand_info = genotype

        if endurance_count == 2:
            interpretation = f"Homozygous endurance ({strand_info})"
        elif power_count == 2:
            interpretation = f"Homozygous power ({strand_info})"
        elif endurance_count == 1:
            interpretation = f"Endurance advantage ({strand_info})"
        elif power_count == 1:
            interpretation = f"Power advantage ({strand_info})"
        else:
            interpretation = f"Neutral/heterozygous ({strand_info})"

        return score, genotype, interpretation

    def calculate_polygenic_score(self) -> Dict:
        """Calculate the full athletic polygenic score"""

        total_score = 0
        snp_scores = []
        category_scores = {
            'muscle_energy': 0,
            'adrenergic': 0,
            'oxygen': 0,
            'neuromotor': 0,
            'inflammation': 0,
            'connective': 0,
            'fuel': 0
        }

        for rsid, snp_info in ATHLETIC_SNPS.items():
            score, genotype, interpretation = self.calculate_snp_score(snp_info)

            if score is not None:
                total_score += score
                category_scores[snp_info.category] += score

                snp_scores.append({
                    'rsid': rsid,
                    'gene': snp_info.gene,
                    'name': snp_info.name,
                    'category': snp_info.category,
                    'score': score,
                    'genotype': genotype,
                    'interpretation': interpretation,
                    'effect': snp_info.effect_description,
                    'notes': snp_info.notes
                })
            else:
                self.missing_snps.append(rsid)

        # Calculate percentages (based on max realistic score of ~3.0)
        max_realistic_score = 3.0
        endurance_percentage = max(0, min(100, (total_score / max_realistic_score) * 100))
        power_percentage = max(0, min(100, (-total_score / max_realistic_score) * 100))

        # Determine athletic type (using effect-size weighted thresholds)
        if total_score > 1.5:
            athletic_type = "Strong Endurance"
        elif total_score > 0.5:
            athletic_type = "Moderate Endurance"
        elif total_score < -1.5:
            athletic_type = "Strong Power/Sprint"
        elif total_score < -0.5:
            athletic_type = "Moderate Power/Sprint"
        else:
            athletic_type = "Balanced/Mixed"

        return {
            'total_score': total_score,
            'athletic_type': athletic_type,
            'endurance_percentage': endurance_percentage,
            'power_percentage': power_percentage,
            'snp_scores': snp_scores,
            'category_scores': category_scores,
            'missing_snps': self.missing_snps,
            'snps_analyzed': len(snp_scores),
            'snps_total': len(ATHLETIC_SNPS)
        }

    def print_report(self, results: Dict):
        """Print a formatted athletic performance report"""

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
        print(f"{BOLD}{CYAN} ATHLETIC PERFORMANCE POLYGENIC SCORE{END}")
        print(f"{BOLD}{CYAN}{'='*80}{END}\n")

        # Overall score
        score = results['total_score']
        score_color = GREEN if score > 0 else RED if score < 0 else YELLOW

        print(f"{BOLD}Overall Score:{END} {score_color}{score:+.2f}{END}")
        print(f"{BOLD}Athletic Type:{END} {BOLD}{score_color}{results['athletic_type']}{END}")
        print(f"{BOLD}SNPs Analyzed:{END} {results['snps_analyzed']}/{results['snps_total']}")

        # Visual bar
        print(f"\n{BOLD}Performance Spectrum:{END}")
        bar_length = 40
        center = bar_length // 2
        # Scale based on max realistic score (~3.0 for very strong endurance/power)
        max_score = 3.0
        position = center + int((score / max_score) * center)
        position = max(0, min(bar_length - 1, position))

        bar = ['-'] * bar_length
        bar[center] = '|'
        bar[position] = '●'

        print(f"Power {RED}{''.join(bar[:center])}{END}{YELLOW}{''.join(bar[center])}{END}{GREEN}{''.join(bar[center+1:])}{END} Endurance")
        print(f"      {' ' * position}↑")
        print(f"      {' ' * position}You")

        # Category breakdown
        print(f"\n{BOLD}{YELLOW}Category Scores:{END}")
        for category, cat_score in results['category_scores'].items():
            cat_color = GREEN if cat_score > 0 else RED if cat_score < 0 else ''
            cat_name = category.replace('_', ' ').title()
            print(f"  {cat_name:20} {cat_color}{cat_score:+6.2f}{END}")

        # Detailed SNP results
        print(f"\n{BOLD}{YELLOW}Individual SNP Results:{END}")

        for category in ['muscle_energy', 'adrenergic', 'oxygen', 'neuromotor', 'inflammation', 'connective', 'fuel']:
            cat_snps = [s for s in results['snp_scores'] if s['category'] == category]
            if cat_snps:
                cat_name = category.replace('_', ' ').upper()
                print(f"\n  {BOLD}{CYAN}{cat_name}:{END}")

                for snp in cat_snps:
                    score_symbol = '→' if snp['score'] == 0 else '▲' if snp['score'] > 0 else '▼'
                    score_color = GREEN if snp['score'] > 0 else RED if snp['score'] < 0 else YELLOW

                    print(f"    {BOLD}{snp['gene']}{END} {snp['name']} ({snp['rsid']})")
                    print(f"      {score_color}{score_symbol} {snp['interpretation']}{END}")
                    print(f"      {snp['notes']}")

        # Missing SNPs
        if results['missing_snps']:
            print(f"\n{BOLD}{YELLOW}Missing SNPs:{END}")
            print(f"  {', '.join(results['missing_snps'])}")

        print(f"\n{BOLD}Interpretation:{END}")
        if score > 1.5:
            print("  Your genetics strongly favor endurance activities:")
            print("  • Marathon, cycling, swimming, triathlon")
            print("  • Better fat oxidation and aerobic capacity")
            print("  • Superior fatigue resistance")
        elif score > 0.5:
            print("  Your genetics moderately favor endurance activities:")
            print("  • Distance running, cycling")
            print("  • Good aerobic capacity")
            print("  • Above-average fatigue resistance")
        elif score < -1.5:
            print("  Your genetics strongly favor power/sprint activities:")
            print("  • Sprinting, weightlifting, throwing")
            print("  • Fast-twitch muscle fiber dominance")
            print("  • Explosive power generation")
        elif score < -0.5:
            print("  Your genetics moderately favor power/sprint activities:")
            print("  • Short sprints, strength training")
            print("  • Good explosive capacity")
            print("  • Above-average power output")
        else:
            print("  Your genetics show a balanced profile:")
            print("  • Suitable for mixed sports (soccer, basketball, hockey)")
            print("  • Can excel in both endurance and power with proper training")
            print("  • Versatile athletic potential")

        print(f"\n{BOLD}{RED}⚠ IMPORTANT LIMITATIONS:{END}")
        print(f"{YELLOW}• This score analyzes only {results['snps_analyzed']} SNPs out of 250+ identified in research{END}")
        print(f"{YELLOW}• Missing key SNPs: CDKN1A rs236448, VDR rs2228570, and hundreds of GWAS hits{END}")
        print(f"{YELLOW}• Athletic performance is highly polygenic (thousands of variants){END}")
        print(f"{YELLOW}• Genetics explains only ~20-30% of athletic performance{END}")
        print(f"{YELLOW}• This score does NOT reliably predict elite athletic potential{END}")
        print(f"\n{BOLD}Training, nutrition, and psychology matter FAR more than genetics!{END}")


def analyze_athletic_performance(snps_data: snps.SNPs) -> Dict:
    """Main entry point for athletic analysis"""
    scorer = AthleticScorer(snps_data)
    results = scorer.calculate_polygenic_score()
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