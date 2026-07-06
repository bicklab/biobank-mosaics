# Population-scale detection of putative embryonic mosaics from biobank CHIP calls

This repository contains the analysis code used to identify putative embryonic mosaic
individuals for canonical clonal hematopoiesis (CHIP) driver genes (e.g. *DNMT3A*, *TET2*)
from large biobank short-read whole-genome sequencing data.

## Rationale

Individuals under 40 years old who carry a high–variant-allele-fraction (VAF ≥ 0.25) somatic
pathogenic variant in a CHIP driver gene are unlikely to have reached that clone size through
age-related clonal hematopoiesis. Such individuals are candidate embryonic mosaics. This code
takes per-cohort CHIP call sets, harmonizes them, classifies putative mosaics by age and VAF,
plots the age–VAF relationship, and estimates prevalence with exact binomial confidence
intervals.

## Inputs

CHIP call sets are generated upstream with the Mutect2 + ANNOVAR pipeline of
Vlasschaert et al., *Blood* 2023 (https://doi.org/10.1182/blood.2022018825), run across
canonical CHIP driver genes. The binomial germline filter is deliberately **omitted** so that
genuine mosaic variants near 50% VAF are retained rather than discarded as germline.

## Method summary

1. Harmonize per-cohort CHIP calls into a common schema (`age, gender, gene, variant, VAF, cohort`).
2. For each gene, classify putative mosaics as `age < 40` and `VAF ≥ 0.25`.
3. Plot age vs. VAF, highlighting putative mosaics.
4. Estimate prevalence per 10,000 people with Clopper–Pearson (exact binomial) confidence intervals.

Candidate mosaics identified this way are then evaluated against alternative explanations
(sequencing artifact, germline/overgrowth-syndrome phenotype, comparison to nearby germline
heterozygous variant VAF, complete blood counts, and co-occurring somatic CHIP mutations).
Those checks rely on cohort-specific clinical tables and are described in the manuscript Methods.

## Data availability

Individual-level biobank data is access-controlled and is **not** included here. The
*All of Us* Research Program and Vanderbilt's BioVU each govern access through their own data
use agreements. Only the analysis code and synthetic example data are shared in this repository.
