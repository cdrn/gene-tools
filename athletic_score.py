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
    notes: Optional[str] = None


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
        notes='T allele = enzyme deficiency, may reduce power output'
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
        notes='A = better endurance, less muscle damage'
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
        notes='G (Gly) = better endurance economy'
    ),

    'rs1042714': AthleticSNP(
        rsid='rs1042714',
        gene='ADRB2',
        name='Gln27Glu',
        category='adrenergic',
        endurance_alleles=['G'],  # Glu
        power_alleles=['C'],      # Gln
        effect_description='Beta-2 receptor downregulation',
        evidence='Bronchodilation, exercise response',
        notes='G (Glu) = better endurance, bronchodilation'
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
        notes='Ancestry-dependent effects, T generally endurance'
    ),

    # ============================================================================
    # INFLAMMATION, RECOVERY & CONNECTIVE TISSUE
    # ============================================================================

    'rs1800795': AthleticSNP(
        rsid='rs1800795',
        gene='IL6',
        name='-174G>C',
        category='inflammation',
        endurance_alleles=['C'],
        power_alleles=['G'],
        effect_description='IL-6 production',
        evidence='Inflammation, recovery',
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
        notes='T = optimal tendon stiffness for endurance'
    ),

    'rs1800012': AthleticSNP(
        rsid='rs1800012',
        gene='COL1A1',
        name='Type I collagen',
        category='connective',
        endurance_alleles=['T'],
        power_alleles=['G'],
        effect_description='Collagen strength',
        evidence='Injury risk, connective tissue',
        notes='T = stronger collagen, endurance advantage'
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
        notes='T = better lactate transport, endurance'
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

    def calculate_snp_score(self, snp_info: AthleticSNP) -> Tuple[Optional[int], str, str]:
        """
        Calculate score for a single SNP
        Returns: (score, genotype, interpretation)
        Score: +2 for homozygous endurance, -2 for homozygous power, 0 for heterozygous
        """
        snp_data = self.snps_data.snps[self.snps_data.snps.index == snp_info.rsid]

        if len(snp_data) == 0 or snp_data['genotype'].isna().all():
            return None, 'Missing', 'Data not available'

        genotype = str(snp_data['genotype'].iloc[0])

        # Count alleles
        endurance_count = 0
        power_count = 0

        for allele in genotype:
            if allele in snp_info.endurance_alleles:
                endurance_count += 1
            elif allele in snp_info.power_alleles:
                power_count += 1

        # Calculate score
        score = endurance_count - power_count

        # Generate interpretation
        if score == 2:
            interpretation = f"Homozygous endurance ({genotype})"
        elif score == -2:
            interpretation = f"Homozygous power ({genotype})"
        elif score == 1:
            interpretation = f"Endurance advantage ({genotype})"
        elif score == -1:
            interpretation = f"Power advantage ({genotype})"
        else:
            interpretation = f"Neutral/heterozygous ({genotype})"

        return score, genotype, interpretation

    def calculate_polygenic_score(self) -> Dict:
        """Calculate the full athletic polygenic score"""

        total_score = 0
        snp_scores = []
        category_scores = {
            'muscle_energy': 0,
            'adrenergic': 0,
            'oxygen': 0,
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

        # Calculate percentages
        max_possible_score = len(ATHLETIC_SNPS) * 2
        endurance_percentage = max(0, (total_score / max_possible_score) * 100)
        power_percentage = max(0, (-total_score / max_possible_score) * 100)

        # Determine athletic type
        if total_score > 5:
            athletic_type = "Strong Endurance"
        elif total_score > 2:
            athletic_type = "Moderate Endurance"
        elif total_score < -5:
            athletic_type = "Strong Power/Sprint"
        elif total_score < -2:
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

        print(f"{BOLD}Overall Score:{END} {score_color}{score:+d}{END}")
        print(f"{BOLD}Athletic Type:{END} {BOLD}{score_color}{results['athletic_type']}{END}")
        print(f"{BOLD}SNPs Analyzed:{END} {results['snps_analyzed']}/{results['snps_total']}")

        # Visual bar
        print(f"\n{BOLD}Performance Spectrum:{END}")
        bar_length = 40
        center = bar_length // 2
        position = center + int((score / (len(ATHLETIC_SNPS) * 2)) * center)
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
            print(f"  {cat_name:20} {cat_color}{cat_score:+3d}{END}")

        # Detailed SNP results
        print(f"\n{BOLD}{YELLOW}Individual SNP Results:{END}")

        for category in ['muscle_energy', 'adrenergic', 'oxygen', 'inflammation', 'connective', 'fuel']:
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
        if score > 5:
            print("  Your genetics strongly favor endurance activities:")
            print("  • Marathon, cycling, swimming, triathlon")
            print("  • Better fat oxidation and aerobic capacity")
            print("  • Superior fatigue resistance")
        elif score > 2:
            print("  Your genetics moderately favor endurance activities:")
            print("  • Distance running, cycling")
            print("  • Good aerobic capacity")
            print("  • Above-average fatigue resistance")
        elif score < -5:
            print("  Your genetics strongly favor power/sprint activities:")
            print("  • Sprinting, weightlifting, throwing")
            print("  • Fast-twitch muscle fiber dominance")
            print("  • Explosive power generation")
        elif score < -2:
            print("  Your genetics moderately favor power/sprint activities:")
            print("  • Short sprints, strength training")
            print("  • Good explosive capacity")
            print("  • Above-average power output")
        else:
            print("  Your genetics show a balanced profile:")
            print("  • Suitable for mixed sports (soccer, basketball)")
            print("  • Can excel in both endurance and power with training")
            print("  • Versatile athletic potential")

        print(f"\n{YELLOW}Note: Genetics is only ~20-30% of athletic performance.{END}")
        print(f"{YELLOW}Training, nutrition, and psychology matter more!{END}")


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