#!/usr/bin/env python3
"""
Alternative entry point for DNA analysis to avoid import issues
"""

import sys
import os

# Add current directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import and run the main analysis
from analyze_dna import analyze_dna, Colors

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
    parser.add_argument('--no-cognitive', action='store_true',
                        help='Skip cognitive/educational attainment analysis')
    parser.add_argument('--no-athletic', action='store_true',
                        help='Skip athletic performance analysis')
    parser.add_argument('--no-haplogroup', action='store_true',
                        help='Skip Y-chromosome haplogroup prediction')
    parser.add_argument('--health-only', action='store_true',
                        help='Only run health-related analyses')
    parser.add_argument('--quick', action='store_true',
                        help='Quick mode: skip polygenic scores (athletic & cognitive)')
    parser.add_argument('--cognitive-threshold', type=float, default=0.001,
                        help='P-value threshold for cognitive score (default: 0.001)')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Show detailed output and statistics')

    args = parser.parse_args()

    # Handle quick mode and health-only mode
    if args.quick:
        args.no_cognitive = True
        args.no_athletic = True
    if args.health_only:
        args.no_cognitive = True
        args.no_athletic = True
        args.no_haplogroup = True

    try:
        # Store args in the analyze_dna module
        import analyze_dna as ad
        ad.CLI_ARGS = args

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