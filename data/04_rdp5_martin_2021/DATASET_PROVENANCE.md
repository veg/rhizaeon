# Dataset Provenance: Martin et al. (2021) Canonical RDP5 Benchmark Suite

## Canonical Citation
> **Martin, D. P., Varsani, A., Roumagnac, P., Botha, G., Maslamoney, S., Schwab, T., Kelz, Z., Kumar, V., & Murrell, B. (2021).**  
> *RDP5: a computer program for analyzing recombination in, and removing signals of recombination from, nucleotide sequence datasets.*  
> **Virus Evolution**, 7(1), veab003.  
> DOI: [10.1093/ve/veab003](https://doi.org/10.1093/ve/veab003)

---

## Dataset Architecture & Original Benchmark Alignments

The canonical benchmark suite evaluated in Martin et al. (2021) Supplementary Table S1 is distributed directly with the RDP5 software package (`RDP5.93Setup.exe`). All alignments are evaluated from the original distribution files without modification:

1. **HIV-1 KAL153 (`Example2(A-J-cons-kal153).fas`):**
   - Retroviral recombinant mosaic across Subtypes A and C in the 5' gag/pol region.
   - Sample size: $N = 9$ taxa.
   - Sequence length: $L = 9,954$ nt.

2. **Potyvirus Genomes (`Example1(PotySeqs).fas`):**
   - Full-length plant potyvirus polyprotein genomes exhibiting multiple recombination crossovers in P1 and NIa-NIb cistrons.
   - Sample size: $N = 25$ taxa.
   - Sequence length: $L = 9,594$ nt.

3. **Foot-and-Mouth Disease Virus (`Example3_FMDV_aligned.fasta`):**
   - Picornaviral field isolates exhibiting inter-serotype capsid mosaicism (VP1/VP2/VP3).
   - Sample size: $N = 91$ taxa.
   - Sequence length: $L = 8,294$ nt.

4. **Pan-Group M HIV-1 (`HIV_Example_aligned.fasta`):**
   - Comprehensive multi-subtype global cohort of Group M HIV-1 complete coding sequences.
   - Sample size: $N = 274$ taxa.
   - Sequence length: $L = 9,557$ nt.

5. **Tomato Yellow Leaf Curl Begomovirus (`tylcv_aligned.fasta`):**
   - Geminivirus circular single-stranded DNA genome with inter-species recombinant isolate AF271234.1 (TYLCV-IS76).
   - Sample size: $N = 6$ taxa.
   - Sequence length: $L = 2,989$ nt.

---

## Published Reference Runtimes (Martin et al. 2021 Table S1)
- **HIV-1 KAL153 ($N=9, L=9{,}954$):** RDP4 = 8.2 s, RDP5 = 1.9 s.
- **Potyvirus Full ($N=25, L=9{,}594$):** RDP4 = 40.2 s, RDP5 = 13.9 s.
- **FMDV Multi-Strain ($N=91, L=8{,}294$):** RDP4 = 15.1 min (907 s), RDP5 = 3.85 min (231 s).
- **Pan-Group M HIV-1 ($N=274, L=9{,}557$):** RDP4 = 4.41 h (15,864 s), RDP5 = 1.02 h (3,662 s).
