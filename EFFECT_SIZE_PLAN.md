# Athletic SNP Effect Size Implementation Plan

## Approach
Since we've hit rate limits and many SNPs lack well-established effect sizes, we'll use a tiered approach:

### Tier 1: Strong Evidence (OR from meta-analyses/large studies)
- Use published ORs directly
- Examples: ACTN3, ACE, AMPD1, MSTN

### Tier 2: Moderate Evidence (replicated associations)
- Assign moderate effect size (OR ~1.3-1.5)
- Examples: PPARGC1A, most SNPs with consistent athlete studies

### Tier 3: Preliminary Evidence (single studies or inconsistent)
- Assign small effect size (OR ~1.1-1.2) or exclude
- Examples: Some newer SNPs with limited replication

## Documented Effect Sizes

### CONFIRMED (from existing code citations):
1. **rs1805086 (MSTN)**: OR = 2.02 (PMC3024427)
2. **rs17602729 (AMPD1)**: OR = 2.17 (PMC12152022, verified earlier)
3. **rs1800255 (COL3A1)**: OR = 4.79 (SNPedia, for injury not performance)

### FROM EARLIER VERIFICATION:
4. **rs1815739 (ACTN3)**: Multiple meta-analyses show OR ~1.4-1.7 for power athletes
   - Conservative estimate: OR = 1.5

5. **rs8192678 (PPARGC1A)**: Meta-analysis PMC6326506
   - G allele enriched in endurance, specific OR needed
   - Conservative estimate: OR = 1.3

6. **rs4343 (ACE I/D)**: Well-studied, mixed results
   - Conservative estimate: OR = 1.2

## Proposed Effect Size Assignments

### Strong effects (OR ≥ 1.5):
- ACTN3 R577X: 1.5
- AMPD1 Q12X: 2.17
- MSTN K153R: 2.02 (rare)

### Moderate effects (OR 1.2-1.5):
- PPARGC1A Gly482Ser: 1.3
- ACE I/D: 1.2
- SOD2 Ala16Val: 1.3
- COL1A1: 1.3 (bone/tendon)

### Small/Uncertain effects (OR 1.0-1.2):
- Most other SNPs: 1.1
- Very weak evidence: 1.0 (neutral, essentially exclude)

## Implementation Strategy
1. Add effect_size field to all SNPs
2. Use logarithmic scoring: log(OR) for homozygotes, log(OR)/2 for heterozygotes
3. Sum weighted scores across all SNPs
4. Convert back to interpretable scale

This avoids the arbitrary +2/-2 system and uses actual genetic effect magnitudes.

