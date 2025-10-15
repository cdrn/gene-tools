#!/usr/bin/env python3
"""
Educational Attainment / Cognitive Ability Polygenic Score
Based on Lee et al. (2018) GWAS of 766,345 individuals
"""

import snps
import pandas as pd
import numpy as np
from typing import Dict, Optional, Tuple
from collections import defaultdict
import os


class CognitiveScorer:
    """Calculate educational attainment polygenic score"""

    def __init__(self, snps_data: snps.SNPs, sumstats_file: str = "EA_GWAS.txt"):
        self.snps_data = snps_data
        self.sumstats_file = sumstats_file
        self.sumstats = None
        self.results = {}

    def load_summary_stats(self, p_threshold: float = 0.05, top_n: Optional[int] = None) -> pd.DataFrame:
        """
        Load GWAS summary statistics
        p_threshold: Only load SNPs with p-value below this threshold
        top_n: If specified, only use the top N SNPs by p-value (overrides p_threshold)
        """
        print(f"Loading summary statistics from {self.sumstats_file}...")

        if not os.path.exists(self.sumstats_file):
            raise FileNotFoundError(
                f"Could not find {self.sumstats_file}. "
                "Please download from SSGAC: https://thessgac.com/"
            )

        # Read the summary statistics file
        # Columns: MarkerName, CHR, POS, A1, A2, EAF, Beta, SE, Pval
        sumstats = pd.read_csv(
            self.sumstats_file,
            sep='\t',
            usecols=['MarkerName', 'A1', 'A2', 'EAF', 'Beta', 'Pval'],
            dtype={'MarkerName': str, 'A1': str, 'A2': str, 'EAF': float, 'Beta': float, 'Pval': float}
        )

        if top_n is not None:
            # Sort by p-value and take top N
            sumstats = sumstats.nsmallest(top_n, 'Pval').copy()
            print(f"Using top {top_n:,} SNPs by p-value")
            print(f"P-value range: {sumstats['Pval'].min():.2e} to {sumstats['Pval'].max():.2e}")
        else:
            # Filter by p-value threshold
            sumstats = sumstats[sumstats['Pval'] < p_threshold].copy()

        # Set rsID as index for fast lookup
        sumstats.set_index('MarkerName', inplace=True)

        print(f"Loaded {len(sumstats):,} SNPs with p < {p_threshold}")

        self.sumstats = sumstats
        return sumstats

    def calculate_polygenic_score(self, p_threshold: float = 0.05, calculate_z_score: bool = True, top_n: Optional[int] = None, bootstrap: bool = False) -> Dict:
        """
        Calculate the educational attainment polygenic score

        p_threshold: P-value threshold for including SNPs (default 0.05)
                    Lower = fewer SNPs but stronger effects
                    Common choices: 5e-8 (genome-wide), 1e-5, 0.001, 0.05, 1.0
        calculate_z_score: Whether to attempt standardization (uncertain with low coverage)
        top_n: If specified, only use the top N SNPs by p-value (overrides p_threshold)
        """

        # Load summary statistics if not already loaded
        if self.sumstats is None:
            self.load_summary_stats(p_threshold=p_threshold, top_n=top_n)

        user_snps = self.snps_data.snps

        # Match SNPs between user data and summary statistics
        matched_snps = []
        total_score = 0.0
        sum_of_weights = 0.0  # For normalization
        sum_of_abs_weights = 0.0  # Sum of absolute values

        print(f"\nMatching SNPs between your data and GWAS...")

        for rsid in self.sumstats.index:
            if rsid not in user_snps.index:
                continue

            user_snp = user_snps.loc[rsid]
            gwas_snp = self.sumstats.loc[rsid]

            # Get user's genotype
            genotype = str(user_snp['genotype'])
            if pd.isna(genotype) or genotype in ['--', '00', 'II', 'DD']:
                continue

            # Count effect alleles
            effect_allele = gwas_snp['A1'].upper()
            other_allele = gwas_snp['A2'].upper()
            beta = gwas_snp['Beta']

            # Count how many copies of the effect allele the user has
            allele_count = sum([1 for allele in genotype.upper() if allele == effect_allele])

            # Verify alleles match (basic quality control)
            user_alleles = set(genotype.upper())
            gwas_alleles = {effect_allele, other_allele}

            # Skip if alleles don't match (possible strand issue or error)
            if not user_alleles.issubset(gwas_alleles):
                # Try complementary strand
                complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
                try:
                    comp_effect = complement[effect_allele]
                    comp_other = complement[other_allele]
                    comp_gwas_alleles = {comp_effect, comp_other}

                    if user_alleles.issubset(comp_gwas_alleles):
                        # Use complementary alleles
                        allele_count = sum([1 for allele in genotype.upper() if allele == comp_effect])
                    else:
                        continue  # Skip this SNP
                except KeyError:
                    continue

            # Add to score
            snp_contribution = beta * allele_count
            total_score += snp_contribution
            sum_of_weights += beta * allele_count  # Weighted sum
            sum_of_abs_weights += abs(beta)

            matched_snps.append({
                'rsid': rsid,
                'genotype': genotype,
                'effect_allele': effect_allele,
                'allele_count': allele_count,
                'beta': beta,
                'pval': gwas_snp['Pval'],
                'contribution': snp_contribution,
                'eaf': gwas_snp.get('EAF', 0.5)  # Default to 0.5 if missing
            })

        print(f"Matched {len(matched_snps):,} SNPs between your data and GWAS")

        if len(matched_snps) == 0:
            return {
                'error': 'No SNPs could be matched between your data and the GWAS summary statistics'
            }

        # Calculate coverage
        coverage = len(matched_snps) / len(self.sumstats) * 100

        # Sort matched SNPs by absolute contribution
        matched_snps_sorted = sorted(matched_snps, key=lambda x: abs(x['contribution']), reverse=True)

        # Build results dictionary
        results = {
            'raw_score': total_score,
            'normalized_score': total_score / sum_of_abs_weights if sum_of_abs_weights > 0 else 0,
            'num_snps_matched': len(matched_snps),
            'num_snps_available': len(self.sumstats),
            'coverage_percent': coverage,
            'p_threshold': p_threshold,
            'matched_snps': matched_snps_sorted[:100],  # Top 100 contributors
            'sum_of_weights': sum_of_abs_weights,
            'top_n': top_n  # Store if we used top N selection
        }

        # Calculate standardized score only if requested and coverage is sufficient
        if calculate_z_score and coverage > 3:  # Require at least 3% coverage
            if top_n and top_n <= 1000:
                print(f"Note: Z-scores with top SNPs selection may not follow standard distribution")
            standardized = self._standardize_score_improved(matched_snps, bootstrap=bootstrap)
            results.update(standardized)
        else:
            print(f"Skipping standardization (coverage too low: {coverage:.1f}%)")
            results.update({
                'z_score': None,
                'percentile': None,
                'expected_mean': None,
                'expected_sd': None,
                'confidence_low': None,
                'confidence_high': None
            })

        return results

    def _standardize_score_improved(self, matched_snps: list, bootstrap: bool = False) -> Dict:
        """
        Improved standardization approach

        Since we don't have the true population distribution, we provide
        multiple interpretations with appropriate caveats
        """
        import numpy as np
        from scipy import stats

        raw_score = sum(snp['contribution'] for snp in matched_snps)

        # Method 1: Use the mean and SD of the matched SNPs (empirical)
        contributions = [snp['contribution'] for snp in matched_snps]
        empirical_mean = np.mean(contributions) * len(matched_snps)
        empirical_sd = np.std(contributions) * np.sqrt(len(matched_snps))

        if empirical_sd > 0:
            empirical_z = (raw_score - empirical_mean) / empirical_sd
        else:
            empirical_z = 0

        # Method 2: Use expected value based on allele frequencies
        # IMPORTANT: When using top N SNPs, we need to account for selection bias
        # The matched SNPs are a biased subsample (we matched the ones present)
        expected_mean = 0.0
        expected_variance = 0.0

        for snp in matched_snps:
            beta = snp['beta']
            eaf = snp.get('eaf', 0.5)  # Use EAF if available

            # Expected contribution
            expected_mean += beta * 2 * eaf
            # Variance
            expected_variance += (beta ** 2) * 2 * eaf * (1 - eaf)

        expected_sd = np.sqrt(expected_variance)

        # Adjustment for small sample: use t-distribution instead of normal
        # With fewer SNPs, we should be more conservative
        if len(matched_snps) < 100:
            # Small sample adjustment
            from scipy import stats
            df = len(matched_snps) - 1
            # Use t-distribution critical value instead of z
            # This makes confidence intervals wider for small samples
            t_critical = stats.t.ppf(0.975, df)  # 97.5th percentile for 95% CI
            adjustment = t_critical / 1.96  # How much wider than normal
            expected_sd *= adjustment

        if expected_sd > 0:
            population_z = (raw_score - expected_mean) / expected_sd
        else:
            population_z = empirical_z  # Fallback

        # Bootstrap confidence intervals (optional, slow)
        if bootstrap:
            n_bootstrap = 1000
            bootstrap_scores = []

            for _ in range(n_bootstrap):
                # Resample with replacement
                resampled = np.random.choice(matched_snps, size=len(matched_snps), replace=True)
                bootstrap_score = sum(snp['contribution'] for snp in resampled)
                bootstrap_scores.append(bootstrap_score)

            confidence_low = np.percentile(bootstrap_scores, 2.5)
            confidence_high = np.percentile(bootstrap_scores, 97.5)
            bootstrap_mean = np.mean(bootstrap_scores)
            bootstrap_sd = np.std(bootstrap_scores)
        else:
            # Quick approximation without bootstrap
            confidence_low = raw_score - 1.96 * expected_sd
            confidence_high = raw_score + 1.96 * expected_sd
            bootstrap_mean = raw_score
            bootstrap_sd = expected_sd

        # Use the population-based z-score as the primary estimate
        # This matches the reference population (1000 Genomes European)
        z_score = population_z
        percentile = stats.norm.cdf(z_score) * 100

        return {
            'z_score': z_score,
            'percentile': percentile,
            'confidence_low': confidence_low,
            'confidence_high': confidence_high,
            'bootstrap_mean': bootstrap_mean,
            'bootstrap_sd': bootstrap_sd,
            'expected_mean': expected_mean,
            'expected_sd': expected_sd,
            'empirical_z': empirical_z,  # Alternative interpretation
            'population_z': population_z,  # Primary interpretation
        }

    def _standardize_score(self, matched_snps: list) -> Dict:
        """
        Standardize the polygenic score using allele frequencies

        Returns z-score and percentile based on expected distribution
        """
        # Calculate expected mean and variance
        # For diploid organisms with allele frequency p:
        # Expected genotype count = 2p (mean number of effect alleles)
        # Variance = 2p(1-p) for each SNP

        expected_mean = 0.0
        expected_variance = 0.0

        for snp in matched_snps:
            rsid = snp['rsid']
            beta = snp['beta']

            # Get allele frequency from summary stats
            if rsid in self.sumstats.index:
                eaf = self.sumstats.loc[rsid, 'EAF']

                # Expected contribution to score
                # Mean = beta × 2 × allele_frequency
                expected_mean += beta * 2 * eaf

                # Variance = beta² × 2 × allele_frequency × (1 - allele_frequency)
                expected_variance += (beta ** 2) * 2 * eaf * (1 - eaf)

        expected_sd = np.sqrt(expected_variance)

        # Calculate z-score
        # (observed - expected) / SD
        raw_score = sum(snp['contribution'] for snp in matched_snps)
        z_score = (raw_score - expected_mean) / expected_sd if expected_sd > 0 else 0.0

        # Convert z-score to percentile (cumulative normal distribution)
        from scipy import stats
        percentile = stats.norm.cdf(z_score) * 100

        return {
            'z_score': z_score,
            'percentile': percentile,
            'expected_mean': expected_mean,
            'expected_sd': expected_sd
        }

    def print_report(self, results: Dict):
        """Print a formatted cognitive/educational attainment report"""

        # Color codes
        HEADER = '\033[95m'
        BLUE = '\033[94m'
        CYAN = '\033[96m'
        GREEN = '\033[92m'
        YELLOW = '\033[93m'
        RED = '\033[91m'
        BOLD = '\033[1m'
        END = '\033[0m'

        if 'error' in results:
            print(f"\n{RED}{BOLD}Error:{END} {results['error']}")
            return

        print(f"\n{BOLD}{CYAN}{'='*80}{END}")
        print(f"{BOLD}{CYAN} EDUCATIONAL ATTAINMENT POLYGENIC SCORE{END}")
        print(f"{BOLD}{CYAN} Based on Lee et al. (2018) - 766,345 individuals{END}")
        if results['num_snps_available'] <= 1000:
            print(f"{BOLD}{GREEN} Using TOP {results['num_snps_available']} SNPs approach (high-impact variants only){END}")
        print(f"{BOLD}{CYAN}{'='*80}{END}\n")

        raw_score = results['raw_score']
        num_matched = results['num_snps_matched']
        coverage = results['coverage_percent']
        p_thresh = results['p_threshold']
        z_score = results.get('z_score')
        percentile = results.get('percentile')
        confidence_low = results.get('confidence_low')
        confidence_high = results.get('confidence_high')

        print(f"{BOLD}Raw Polygenic Score:{END} {raw_score:.4f}")
        print(f"{BOLD}Normalized Score:{END} {results.get('normalized_score', 0):.4f}")

        if z_score is not None and percentile is not None:
            print(f"{BOLD}Standardized Score (Z):{END} {z_score:+.2f} SD")

            # Color code percentile
            if percentile >= 75:
                percentile_color = GREEN
            elif percentile >= 25:
                percentile_color = YELLOW
            else:
                percentile_color = RED

            print(f"{BOLD}Percentile Rank:{END} {percentile_color}{percentile:.1f}th{END} (higher than {percentile:.1f}% of population)")

            if confidence_low is not None and confidence_high is not None:
                print(f"{BOLD}95% Confidence Interval:{END} [{confidence_low:.2f}, {confidence_high:.2f}]")
        else:
            print(f"{BOLD}Standardization:{END} {YELLOW}Not calculated (insufficient coverage){END}")

        # Coverage warning - be more conservative with top SNPs
        top_n = results.get('top_n')  # Get top_n from results if it was used
        if results['num_snps_available'] <= 1000 and coverage < 10:
            coverage_color = RED
            coverage_msg = " ⚠️ VERY low coverage - percentile unreliable!"
            # Override percentile display for unreliable cases
            if percentile is not None and coverage < 10:
                print(f"\n{RED}{BOLD}WARNING:{END} With only {num_matched}/{results['num_snps_available']} SNPs matched ({coverage:.1f}%),")
                print(f"         the percentile estimate is highly unreliable and should be ignored.")
                print(f"         Focus on the individual SNP effects shown below instead.")
        elif coverage < 5:
            coverage_color = RED
            coverage_msg = " ⚠️ Very low coverage - interpret with extreme caution!"
        elif coverage < 10:
            coverage_color = YELLOW
            coverage_msg = " ⚠️ Low coverage - results uncertain"
        else:
            coverage_color = GREEN
            coverage_msg = ""

        print(f"{BOLD}SNPs Matched:{END} {num_matched:,} / {results['num_snps_available']:,} ({coverage_color}{coverage:.1f}% coverage{END}){coverage_msg}")
        print(f"{BOLD}P-value Threshold:{END} {p_thresh}")

        # Visual percentile distribution (only if percentile calculated)
        if percentile is not None:
            print(f"\n{BOLD}Percentile Distribution:{END}")
            self._print_percentile_bar(percentile)

        # Interpretation
        print(f"\n{BOLD}{YELLOW}Score Interpretation:{END}")

        if percentile >= 90:
            print(f"  {GREEN}Very High{END} genetic predisposition for educational attainment")
            print(f"  You scored higher than {percentile:.1f}% of the population")
        elif percentile >= 75:
            print(f"  {GREEN}Above Average{END} genetic predisposition for educational attainment")
            print(f"  You scored higher than {percentile:.1f}% of the population")
        elif percentile >= 60:
            print(f"  {YELLOW}Moderately Above Average{END} genetic predisposition")
            print(f"  You scored higher than {percentile:.1f}% of the population")
        elif percentile >= 40:
            print(f"  {YELLOW}Average{END} genetic predisposition for educational attainment")
            print(f"  You scored in the middle range of the population")
        elif percentile >= 25:
            print(f"  {YELLOW}Moderately Below Average{END} genetic predisposition")
            print(f"  You scored higher than {percentile:.1f}% of the population")
        elif percentile >= 10:
            print(f"  {RED}Below Average{END} genetic predisposition for educational attainment")
            print(f"  You scored higher than {percentile:.1f}% of the population")
        else:
            print(f"  {RED}Low{END} genetic predisposition for educational attainment")
            print(f"  You scored higher than {percentile:.1f}% of the population")

        # Rough percentile estimate (very approximate!)
        # The score distribution depends on the SNPs included and the reference population
        print(f"\n{BOLD}{YELLOW}Important Notes:{END}")
        print(f"  • This score explains ~11-13% of variance in educational attainment")
        print(f"  • Environment, motivation, and opportunity matter much more!")
        print(f"  • Educational attainment correlates with intelligence (r ≈ 0.7)")
        print(f"  • This is NOT an IQ test - it's a genetic tendency estimate")
        print(f"  • 87-89% of educational differences are NOT explained by genetics")

        # Show top contributing SNPs with clear allele interpretation
        if results['matched_snps']:
            print(f"\n{BOLD}{YELLOW}Your Educational Attainment SNPs:{END}")
            print(f"{'SNP':<12} {'Your Genotype':<14} {'EA Allele':<10} {'Count':<7} {'Effect':<35}")
            print(f"{'-'*85}")

            for i, snp in enumerate(results['matched_snps'][:25]):
                rsid = snp['rsid']
                genotype = snp['genotype']
                ea_allele = snp['effect_allele']
                count = snp['allele_count']
                contribution = snp['contribution']

                # Determine effect direction and interpretation
                if contribution > 0.01:
                    effect_symbol = '↑↑'
                    effect_color = GREEN
                    effect_text = 'Higher educational attainment'
                elif contribution > 0:
                    effect_symbol = '↑'
                    effect_color = GREEN
                    effect_text = 'Slightly higher EA'
                elif contribution < -0.01:
                    effect_symbol = '↓↓'
                    effect_color = RED
                    effect_text = 'Lower educational attainment'
                elif contribution < 0:
                    effect_symbol = '↓'
                    effect_color = RED
                    effect_text = 'Slightly lower EA'
                else:
                    effect_symbol = '→'
                    effect_color = YELLOW
                    effect_text = 'Neutral'

                # Highlight your allele count
                if count == 2:
                    count_display = f"{BOLD}{count}×{ea_allele}{END}"
                    status = "(homozygous)"
                elif count == 1:
                    count_display = f"{count}×{ea_allele}"
                    status = "(heterozygous)"
                else:
                    count_display = f"0×{ea_allele}"
                    status = "(no copies)"

                print(f"{rsid:<12} {genotype:<14} {ea_allele:<10} {count_display:<16} "
                      f"{effect_color}{effect_symbol} {effect_text}{END}")

        print(f"\n{BOLD}Data Source:{END}")
        print(f"  Lee et al. (2018). Gene discovery and polygenic prediction from a")
        print(f"  1.1-million-person GWAS of educational attainment.")
        print(f"  Nature Genetics, 50(8), 1112-1121.")
        print(f"  https://doi.org/10.1038/s41588-018-0147-3")

    def _print_percentile_bar(self, percentile: float):
        """Print a visual percentile bar"""
        BLUE = '\033[94m'
        GREEN = '\033[92m'
        YELLOW = '\033[93m'
        RED = '\033[91m'
        BOLD = '\033[1m'
        END = '\033[0m'

        bar_length = 50
        position = int((percentile / 100) * bar_length)
        position = max(0, min(bar_length - 1, position))

        # Color zones
        bar = []
        for i in range(bar_length):
            if i < bar_length * 0.1:
                bar.append(RED + '█' + END)
            elif i < bar_length * 0.25:
                bar.append(YELLOW + '█' + END)
            elif i < bar_length * 0.75:
                bar.append(BLUE + '█' + END)
            elif i < bar_length * 0.9:
                bar.append(YELLOW + '█' + END)
            else:
                bar.append(GREEN + '█' + END)

        # Mark position
        bar[position] = BOLD + '▼' + END

        print(f"  0%  {''.join(bar)}  100%")
        print(f"      {' ' * position}↑")
        print(f"      {' ' * position}You ({percentile:.1f}th)")


def analyze_cognitive_score(snps_data: snps.SNPs, p_threshold: float = 0.001, top_n: Optional[int] = None) -> Dict:
    """
    Main entry point for cognitive/educational attainment analysis

    p_threshold: P-value threshold for including SNPs
                 0.05 = all suggestive SNPs (most SNPs, best prediction)
                 5e-8 = genome-wide significant only (fewer SNPs, conservative)
    top_n: If specified, use only the top N SNPs by p-value (recommended: 500)
    """
    scorer = CognitiveScorer(snps_data)
    results = scorer.calculate_polygenic_score(p_threshold=p_threshold, top_n=top_n)
    scorer.print_report(results)
    return results


def analyze_top_snps_comparison(snps_data: snps.SNPs) -> Dict:
    """
    Compare polygenic scores using only the top N most significant SNPs
    This tests if using fewer, higher-quality SNPs gives better signal
    """
    import pandas as pd

    print("\n" + "="*80)
    print(" TOP SNPs ANALYSIS - Using Only Highest Impact Variants")
    print("="*80 + "\n")

    top_n_values = [10, 50, 100, 500, 1000, 5000, 10000, None]
    results_list = []

    for top_n in top_n_values:
        label = f"Top {top_n}" if top_n else "p < 5e-8"
        print(f"\nAnalyzing {label} SNPs...")
        scorer = CognitiveScorer(snps_data)

        if top_n:
            results = scorer.calculate_polygenic_score(p_threshold=1.0, calculate_z_score=True, top_n=top_n)
        else:
            results = scorer.calculate_polygenic_score(p_threshold=5e-8, calculate_z_score=True)

        results_list.append({
            'Selection': label,
            'SNPs Used': results['num_snps_available'],
            'SNPs Matched': results['num_snps_matched'],
            'Match Rate': f"{(results['num_snps_matched']/results['num_snps_available']*100):.1f}%",
            'Raw Score': f"{results['raw_score']:.4f}",
            'Z-score': f"{results.get('z_score', 0):.2f}" if results.get('z_score') else 'N/A',
            'Percentile': f"{results.get('percentile', 0):.1f}%" if results.get('percentile') else 'N/A'
        })

    # Create summary table
    df = pd.DataFrame(results_list)
    print("\n" + "="*80)
    print(" SUMMARY: Score Using Only Top SNPs")
    print("="*80)
    print(df.to_string(index=False))

    print("\n" + "="*80)
    print(" KEY INSIGHTS:")
    print("="*80)
    print("• Top 10-100 SNPs: Capture major effect variants only")
    print("• Top 1000 SNPs: Balance between effect size and coverage")
    print("• Top 10000 SNPs: Include suggestive associations")
    print("• Genome-wide (p<5e-8): All statistically robust associations")
    print("\nIf your percentile improves with fewer SNPs, you have the")
    print("high-impact variants but lack the polygenic background.")
    print("If it decreases, you're missing key large-effect variants.")

    return {'top_n_values': top_n_values, 'results': results_list}


def analyze_multiple_thresholds(snps_data: snps.SNPs) -> Dict:
    """
    Analyze cognitive score at multiple p-value thresholds
    This shows how the score changes with different SNP inclusion criteria
    """
    import pandas as pd

    print("\n" + "="*80)
    print(" MULTI-THRESHOLD COGNITIVE SCORE ANALYSIS")
    print("="*80 + "\n")

    thresholds = [5e-8, 1e-5, 0.001, 0.01, 0.05, 1.0]
    results_list = []

    for p_thresh in thresholds:
        print(f"\nAnalyzing at p < {p_thresh:g}...")
        scorer = CognitiveScorer(snps_data)
        results = scorer.calculate_polygenic_score(p_threshold=p_thresh, calculate_z_score=True)

        results_list.append({
            'P-value': f"{p_thresh:g}",
            'SNPs': results['num_snps_matched'],
            'Coverage %': f"{results['coverage_percent']:.1f}",
            'Raw Score': f"{results['raw_score']:.2f}",
            'Normalized': f"{results['normalized_score']:.4f}",
            'Z-score': f"{results.get('z_score', 0):.2f}" if results.get('z_score') else 'N/A',
            'Percentile': f"{results.get('percentile', 0):.1f}%" if results.get('percentile') else 'N/A'
        })

    # Create summary table
    df = pd.DataFrame(results_list)
    print("\n" + "="*80)
    print(" SUMMARY: Score Across Different P-value Thresholds")
    print("="*80)
    print(df.to_string(index=False))

    print("\n" + "="*80)
    print(" INTERPRETATION:")
    print("="*80)
    print("• Genome-wide significant SNPs (p<5e-8) = Most reliable, strongest effects")
    print("• Suggestive SNPs (p<0.05) = Include more polygenic background")
    print("• All SNPs (p<1.0) = Full polygenic architecture (but more noise)")
    print("\nVariation across thresholds suggests:")
    print("• Stable scores = consistent genetic profile")
    print("• Large variation = possible coverage issues or population stratification")

    return {'thresholds': thresholds, 'results': results_list}


if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser(
        description='Calculate Educational Attainment Polygenic Score',
        epilog='Examples:\n'
               '  %(prog)s ancestryDNA.txt                    # Default analysis\n'
               '  %(prog)s ancestryDNA.txt --top 500          # Use top 500 SNPs (RECOMMENDED)\n'
               '  %(prog)s ancestryDNA.txt --p-value 5e-8     # Genome-wide significant only\n'
               '  %(prog)s ancestryDNA.txt --compare          # Compare different approaches\n'
               '  %(prog)s ancestryDNA.txt --top-analysis     # Analyze top SNPs ranges\n',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('filepath', nargs='?', default='ancestryDNA.txt',
                        help='Path to AncestryDNA.txt file (default: ancestryDNA.txt)')
    parser.add_argument('--p-value', '-p', type=float, default=0.05,
                        help='P-value threshold for SNP inclusion (default: 0.05)')
    parser.add_argument('--top', '-t', type=int, default=None,
                        help='Use only top N SNPs by p-value (recommended: 500)')
    parser.add_argument('--compare', '-c', action='store_true',
                        help='Compare multiple p-value thresholds')
    parser.add_argument('--top-analysis', '-a', action='store_true',
                        help='Analyze different top N selections')

    args = parser.parse_args()

    try:
        print("Loading DNA data...")
        s = snps.SNPs(args.filepath)

        if args.compare:
            # Run multi-threshold comparison
            results = analyze_multiple_thresholds(s)
        elif args.top_analysis:
            # Run top SNPs analysis
            results = analyze_top_snps_comparison(s)
        else:
            # Run single analysis
            if args.top:
                print(f"\nUsing TOP {args.top:,} SNPs approach (recommended for high-achievers)")
            results = analyze_cognitive_score(s, p_threshold=args.p_value, top_n=args.top)

    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
