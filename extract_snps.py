#!/usr/bin/env python3
"""Extract all SNP definitions from analyze_dna.py for verification"""

import re

def extract_snp_data(filepath):
    with open(filepath, 'r') as f:
        content = f.read()

    # Split by SNP entries
    snp_blocks = re.split(r"(\s+'rs\d+(?:_\w+)?'\s*:\s*\{)", content)

    snps = []
    for i in range(1, len(snp_blocks), 2):
        if i+1 < len(snp_blocks):
            rsid_match = re.search(r"'(rs\d+(?:_\w+)?)'", snp_blocks[i])
            if rsid_match:
                rsid = rsid_match.group(1)
                # Get the block content until the next SNP or end
                block_content = snp_blocks[i+1].split("'rs")[0]

                # Extract gene
                gene_match = re.search(r"'gene':\s*'([^']+)'", block_content)
                gene = gene_match.group(1) if gene_match else "Unknown"

                # Extract all genotype interpretations
                genotypes = {}
                for gt in ['CC', 'CT', 'TT', 'GG', 'AG', 'AA', 'CG', 'GC', 'AC', 'CA', 'GT', 'TG', 'GA', 'TC']:
                    pattern = rf"'{gt}':\s*'([^']+)'"
                    match = re.search(pattern, block_content)
                    if match:
                        genotypes[gt] = match.group(1)

                if genotypes:  # Only include if we found genotypes
                    snps.append({
                        'rsid': rsid,
                        'gene': gene,
                        'genotypes': genotypes
                    })

    return snps

if __name__ == "__main__":
    snps = extract_snp_data('/Users/cdrn/Code/gene-analysis/analyze_dna.py')

    print(f"Total SNPs with genotype data: {len(snps)}\n")

    for snp in snps:
        print(f"{snp['rsid']} ({snp['gene']})")
        for gt, interpretation in sorted(snp['genotypes'].items()):
            print(f"  {gt}: {interpretation[:80]}{'...' if len(interpretation) > 80 else ''}")
        print()
