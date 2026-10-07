# Changelog

All notable changes to the GeneMolX WGS pipeline are recorded here.
Versions follow [Semantic Versioning](https://semver.org/).

## [0.2.0] — 2026-10-06

### Added
- ACMG 2021 carrier screening panel (tier 3: 97 autosomal recessive + 16 X-linked genes,
  carrier frequency >=1/200) as `resources/acmg/ACMG_carrier_2021_tier3_genes.tsv`. Carrier
  findings are now reported against this panel and separated from carriers in recessive genes
  outside it, which no guideline recommends screening.
- `PGT`, `PID` and `PS` columns in the curated tables, carrying HaplotypeCaller's physical
  phasing. These distinguish two variants in cis from in trans — the difference between a
  benign haplotype and a compound-heterozygous genotype.
- `ACMG_Carrier_Panel` and `Carrier_Condition` columns in the tiered clinical output.
- Limitation stating that SMN1 copy number and FMR1 repeat expansion cannot be assessed from
  short-read WGS, so their absence from the results is not a negative screening result.

### Changed
- Curated tables: 48 -> 51 columns. Tiered clinical output: 51 -> 53 columns. All additions are
  appended; no existing column moved, was renamed, or changed type.
- Gene inheritance is resolved in the order ACMG SF v3.3 -> ACMG carrier panel -> GenCC ->
  Funcotator.

## [0.1.0] — 2026-09-16

First release. Germline whole-genome short-variant analysis, FASTQ to clinical report.
Validated end to end: full chain executed on a combined test dataset, and steps 4-6 on a
real 30x genome. See the Roadmap in the README for what is planned beyond this release.

### Added
- `--version` / `-V` on every script, reporting pipeline version 0.1.0.
- Sample acceptance gate in `GMX_clinical_variants.sh`: mean coverage ≥ 20x, ≥ 80% of the
  genome at ≥ 20x, contamination ≤ 3%. A failing sample is reported with a prominent banner
  and no variant is classified.
- `DP_outlier` quality flag for variants whose depth far exceeds the genome mean — the
  signature of a collapsed repeat.
- `GMX_fastqc.sh` — first pipeline step. Validates raw FASTQ before anything else runs
  (record structure, sequence/quality length agreement, quality encoding, optional
  decompression integrity) and exits non-zero on bad input. Accepts provider deliveries as
  they come: recursive directories, `.fastq.gz` / `.fq.gz` / `.bz2` / uncompressed, and the
  `_R1`/`_1`/`.R1.` pair conventions. Reports detected pairing and platform, runs FastQC
  with `--nogroup`, aggregates with MultiQC (including Picard and VerifyBamID2 metrics via
  `--include`), prints WARN/FAIL modules, and derives guidance for the trimming parameters.
- `GMX_trim_reads.sh` — Trimmomatic PE adapter and quality trimming.
- `GMX_variant_calling.sh` — BWA-MEM → MarkDuplicatesSpark → BQSR → HaplotypeCaller (GVCF) →
  GenotypeGVCFs, with VerifyBamID2 contamination passed to the caller.
- `GMX_annotation.sh` — GATK hard filters, normalization (multi-allelic split, left-align,
  `*` allele removal), Funcotator with MANE-preferred transcripts, and a 48-column curated
  clinical table including the GATK quality annotations.
- `GMX_clinical_variants.sh` — evidence-tiered clinical output, gnomAD v4.1 filtering allele
  frequency by remote lookup, automated ACMG/AMP classification (ClinGen point system), and
  text + HTML reports.
- `GMX_download_references.sh` — reference genome, known sites, Funcotator bundle, ClinVar,
  AlphaMissense, REVEL, SpliceAI, GenCC, gnomAD constraint, transcript exon structure.
- `environment.yml` — pinned conda environment (`genemolx-wgs`).
- `resources/` — ACMG SF v3.3 gene tables, adapter sequences, GenCC gene-inheritance table,
  gnomAD gene-constraint table.
