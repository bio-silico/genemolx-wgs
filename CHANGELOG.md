# Changelog

All notable changes to the GeneMolX WGS pipeline are recorded here.
Versions follow [Semantic Versioning](https://semver.org/).

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
