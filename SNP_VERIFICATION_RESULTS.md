# SNP Verification Results

## Status: IN PROGRESS
Checking all 166 SNPs across analyze_dna.py, all 37 SNPs in athletic_score.py

---

## ERRORS FOUND

### 1. ❌ rs1800012 (COL1A1 Sp1) - WRONG in analyze_dna.py line 1550-1556

**Current (INCORRECT):**
- GG: Normal collagen
- GT: Stronger collagen (Sp1 site)
- TT: Strongest collagen, better bone density

**Should be:**
- **GG: Normal/best** (highest BMD, lowest fracture risk)
- **GT: Intermediate** (moderately reduced BMD)
- **TT: Worst** (2.97x vertebral fracture risk, lowest BMD, imbalanced collagen)

**Evidence:**
- T allele causes 3-fold higher transcription, creating 2.3:1 instead of 2:1 alpha-1:alpha-2 ratio
- TT genotype associated with reduced BMD and increased fracture risk (OR=1.78-2.97)
- Sources: [PMC5650439](https://pmc.ncbi.nlm.nih.gov/articles/PMC5650439/), [Nature Genetics](https://www.nature.com/articles/ng1096-203)

---

## VERIFIED CORRECT ✅

### Pharmacogenomics (all checked)
1. ✅ rs1799853 (CYP2C9*2) - CORRECT
2. ✅ rs1057910 (CYP2C9*3) - CORRECT
3. ✅ rs4244285 (CYP2C19*2) - CORRECT
4. ✅ rs4986893 (CYP2C19*3) - CORRECT
5. ✅ rs12248560 (CYP2C19*17) - CORRECT
6. ✅ rs1065852 (CYP2D6*4) - CORRECT
7. ✅ rs3745274 (CYP2B6*6) - CORRECT
8. ✅ rs1801133 (MTHFR C677T) - CORRECT
9. ✅ rs1801131 (MTHFR A1298C) - CORRECT

### Disease Risk SNPs (verified)
10. ✅ rs429358 (APOE ε4) - CORRECT (C=risk)
11. ✅ rs7412 (APOE ε2) - CORRECT (T=protective)
12. ✅ rs1333049 (9p21 CAD) - CORRECT (C=risk, G=protective)
13. ✅ rs6025 (Factor V Leiden) - CORRECT (A=risk, G=normal)
14. ✅ rs1799963 (Prothrombin) - CORRECT (A=risk, G=normal)
15. ✅ rs1800562 (HFE C282Y) - CORRECT (A=risk, G=normal)
16. ✅ rs7903146 (TCF7L2) - CORRECT (T=risk, C=normal)
17. ✅ rs738409 (PNPLA3) - CORRECT (G=risk, C=normal)
18. ✅ rs662799 (APOA5) - CORRECT (G=risk, A=protective)
19. ✅ rs3764261 (CETP) - CORRECT (A=beneficial, C=typical)
20. ✅ rs9536314 (KLOTHO) - CORRECT
21. ✅ rs4293393 (UMOD) - CORRECT (A=risk on complement strand)

### Athletic Performance SNPs (athletic_score.py) - ALL VERIFIED ✅
22. ✅ rs1815739 (ACTN3) - CORRECT
23. ✅ rs4343 (ACE) - CORRECT
24. ✅ rs8192678 (PPARGC1A) - CORRECT
25. ✅ rs17602729 (AMPD1) - CORRECT
26. ✅ rs4880 (SOD2) - CORRECT
27. ✅ rs1800012 (COL1A1) - CORRECT in athletic_score.py
28. ✅ rs35767 (IGF1) - CORRECT
29. ✅ rs6265 (BDNF) - CORRECT
30. ✅ rs2228570 (VDR FokI) - CORRECT
31-58. ✅ All remaining 28 athletic SNPs - CORRECT

---

## CHECKING IN PROGRESS

Currently verifying remaining 137 SNPs in analyze_dna.py...

