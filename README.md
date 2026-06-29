**Variant Priositisation Pipeline**
This is a nextflow workflow which takes variant called information in the form of a .VCF file and outputs a list of scored variants as a .tsv file in the results/ folder.
This workflow was developed using 1000 Genome Project chromome 22 data which can be found here https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/. This data is not linked to phenotypes due to data security protocols, thus scoring is only related to the characteristis given to that variant from snpeff, ClinVar and gnomAD annotations and is not associated with a phenotype. The workflow was also developed on a limited computer system so RAM and cpu constraints have been built into the nextflow.config file.
The pipeline is as follows:
```
Multi-sample VCF
      │
      ▼
1. extract_sample   — bcftools view: per-sample VCF extraction
      │
      ▼
2. filter_variants  — Hard quality filters many missing values in the 1000 genome data so only AF and GQ used to filter as place holders.
      │
      ▼
3. snpeff           — Functional annotation (GRCh37.75)
      │
      ▼
4. clinvar          — ClinVar annotation (CLNSIG, CLNDN)
      │
      ▼
5. gnomad           — gnomAD AF annotation + rare variant filter (AF < 0.01)
      │
      ▼
6. make_tsv         — bcftools query to TSV format
      │
      ▼
7. score            — ACMG/AMP-informed Python scoring + tiering
      │
      ▼
results/ranked/
  ├── {sample}_ranked.tsv           # all variants with scores
  └── {sample}_ranked.priority.tsv # Tier 1 + 2 only
```

## Filtering Strategy
 
### Quality filters (filter_variants)
GQ and MQ missing from the test data however ideal scenario:
| FORMAT/GQ | > 20 | Genotype quality
| FORMAT/DP | > 10 | Read depth
| MQ | > 40 | Mapping quality
| INFO/AF | < 1 | (AF=1 present in all samples)
 
### Population frequency filter (gnomad)
Variants with gnomAD AF ≥ 0.001 (0.1%) are removed prior to scoring. This threshold aligns with PM2 evidence in the ACMG/AMP framework — variants above this frequency are unlikely to be causative for rare Mendelian disease. Variants absent from gnomAD are retained and scored positively.

## Scoring
 
### ACMG/AMP criteria implemented
 
The scoring script (`scripts/score_variants.py`) evaluates each variant across four dimensions, each mapping to one or more ACMG/AMP criteria:
 
| Dimension | Criteria | Data source |
|-----------|----------|-------------|
| Population frequency | PM2, BA1, BS1 | gnomAD AF |
| Functional impact | PVS1 (partial), PM4 | SnpEff IMPACT |
| Molecular consequence | PVS1 (partial), PM4, BP7 | SnpEff effect term |
| Clinical significance | PS1, PP5, BP6 | ClinVar CLNSIG |
 
Each criterion contributes a score and a reason string, making every prioritisation decision visible.

### Criteria NOT implemented (require data beyond VCF)
 
The following ACMG/AMP criteria cannot be evaluated from VCF data alone and are not implemented:
 
- **PS2, PM6** — require confirmed de novo (trio sequencing)
- **PS3, BS3** — require functional study results
- **PS4** — requires case-control prevalence data
- **PM1** — requires protein domain / hotspot databases
- **PM3, BP2** — require phase and family data
- **PM5** — requires amino-acid-level ClinVar lookup
- **PP1, BS4** — require family segregation data
- **PP2, BP1** — require gene-level constraint data
- **PP4, BP5** — require phenotype / HPO terms
- **BS2** — requires inheritance mode and carrier data
- **BP3** — requires repeat region annotation

### Score tiers
Tier 1 — ClinVar Confirmed | ClinVar Pathogenic exact match | Urgent clinical review |
Tier 1 — High Priority | Total score ≥ 15 | Clinical review |
Tier 2 — Moderate Priority | Total score 8–14 | Clinical review |
Tier 3 — Low Priority | Total score 3–7 | Deprioritised |
Tier 4 — Likely Benign | Total score < 3 | Not reported |


### Reference files required
Multi-sample VCF + TBI User-provided (e.g. 1000 Genomes chr22)
ClinVar VCF + TBI [NCBI ClinVar](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/)
gnomAD exomes VCF + TBI [gnomAD v4.1](https://gnomad.broadinstitute.org/downloads)

## Installation
 
```bash
git clone https://github.com/YOUR_USERNAME/variant-prioritisation-pipeline.git
cd variant-prioritisation-pipeline
 
# Create conda environment
conda env create -f envs/environment.yml
conda activate variantDx
 
# Download SnpEff database
snpEff download GRCh37.75
```

## Usage
 
```bash
nextflow run main.nf \
    --vcf     /path/to/cohort.vcf.gz \
    --clinvar /path/to/clinvar.vcf.gz \
    --gnomad  /path/to/gnomad.exomes.vcf.bgz \
    --samples test_data/test_samples.txt \
    --outdir  results
```


### Test with 1000 Genomes chr22 data
The pipeline was developed and tested using:
- 1000 Genomes Project chr22 VCF (public, non-clinical data)
- Five unrelated samples: HG00098, HG00100, HG00106, HG00112, HG00114
These are healthy population controls — no clinical phenotype is associated with these samples and variant scores should not be interpreted as clinically meaningful.
 
### `{sample}_ranked.tsv`
Full scored variant table, sorted by total score descending.
 
Column | Description
chr, pos, ref, alt | Variant coordinates 
gene | Gene symbol (SnpEff) 
hgvs_c, hgvs_p | HGVS nomenclature 
effect, impact | SnpEff consequence terms 
af | gnomAD allele frequency 
clin_sig, disease | ClinVar classification and disease name 
af_score, impact_score, effect_score, clinvar_score | Per-dimension scores 
total_score | Sum of all dimension scores 
tier | Prioritisation tier 
*_reason | Verbose string for each scoring decision
 

## Limitations
 
- Supports **GRCh37 only** — reference files and SnpEff database are hg19/GRCh37
- **Single-sample scoring** — variants are scored independently with no family context
- **No phenotype integration** — all variants are scored equally regardless of gene-disease relevance; a gene panel BED file (`--panel`) would be required for phenotype-driven prioritisation
- **Partial ACMG/AMP implementation** — ~7–9 of 28 criteria are computable from VCF data; this is an inherent limitation of computational classification
- **Not validated for clinical use** — thresholds and weights are based on the ACMG/AMP framework but have not been prospectively validated against a clinical variant dataset
---
