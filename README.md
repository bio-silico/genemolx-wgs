# GeneMolX WGS Pipeline

Your whole-genome sequencing data, from raw FASTQ to a clinical variant report you can read.

Built on GATK4 best practices. Interpretation follows ACMG/AMP 2015 and current ClinGen
guidance. Works with FASTQ from any provider — Illumina, MGI/DNBSEQ, gzip, bzip2, one lane
or many.

> **Whole-genome sequencing only.** This pipeline is designed for WGS. Exome and targeted
> panels are not supported: coverage metrics are computed genome-wide and variant calling runs
> without target intervals, so both would be wrong on captured data.

> **Research and informational use only. Not a diagnostic device.**
> Nothing here is a clinical result. Any finding you intend to act on must be confirmed by an
> accredited laboratory and interpreted by a qualified clinician.

---

## What you get

| Output | What it is |
|---|---|
| `clinical_report.html` | the report — open it in a browser, print it |
| `clinical_report.txt` | same report, plain text |
| `clinical_tiers.tsv` | every candidate variant, 51 columns, one row each, with the evidence behind it |
| `curated_snps_clinical.tsv` / `curated_indels_clinical.tsv` | **all** your variants (~4.5 M), 48 columns, for your own filtering |
| `multiqc_report.html` | sequencing quality, all metrics in one page |

The report answers, in this order: is the data good enough · anything in the 84 ACMG
secondary-findings genes · known pathogenic variants · carrier status · rare candidates
worth a look · what needs human review · what this analysis cannot see.

---

## Requirements

| | |
|---|---|
| OS | Linux, bash |
| CPU | 8 cores minimum, 16+ recommended |
| RAM | 32 GB minimum, 64 GB recommended |
| Disk | **~320 GB** for one 30x genome: 80 GB reference data + 240 GB working files |
| Software | conda or mamba |
| Network | needed for reference downloads, and during annotation (gnomAD 2.1 is read remotely) |

Time for one 30x genome, on 24 cores: **about 24 hours**, most of it variant calling.

---

## 1. Install

```bash
git clone https://github.com/bio-silico/genemolx-wgs.git
cd genemolx-wgs
conda env create -f environment.yml
conda activate genemolx-wgs
```

Every script activates the environment itself, so you can also run them from cron or a
scheduler without activating first.

## 2. Get the reference data

```bash
./GMX_download_references.sh /path/to/references --all
```

One command, then wait. Roughly 80 GB and several hours on a normal connection.

| Section | Size | Flag |
|---|---|---|
| GRCh38 genome + indexes | 20 GB | `--reference` |
| BQSR known sites (dbSNP, Mills, known indels) | included above | `--known-sites` |
| Funcotator bundle (Gencode v43, gnomAD, …) | 12 GB | `--funcotator` |
| ClinVar (current release) | 1.2 GB | `--clinvar` |
| AlphaMissense | 410 MB | `--alphamissense` |
| REVEL | 335 MB | `--revel` |
| SpliceAI (Ensembl MANE, SNVs) | 8.3 GB | `--spliceai` |
| GenCC gene-disease validity | 26 MB | `--gencc` |
| gnomAD gene constraint | 91 MB | `--constraint` |
| Transcript exon structure | built locally | `--exon-structure` |

Resume-safe: rerun the same command and finished downloads are skipped.

## 3. Run it

> **Before you run:** every input is passed as an argument — nothing is hard-coded to a
> particular sample. Set the paths below to your own reference directory, read files and
> output directory, and use your own sample name.

```bash
REF=/path/to/references/hg38
FUNC=/path/to/references/funcotator_dataSources.v1.8.hg38.20230908g
OUT=/path/to/mysample

# 1. Check the raw data — do this first, it validates your files and tells you how to trim
./GMX_fastqc.sh $OUT/qc_raw /path/to/raw_reads --threads 8

# 2. Trim
./GMX_trim_reads.sh raw_R1.fq.gz raw_R2.fq.gz $OUT/reads MYSAMPLE '' 8

# 3. Check the trimmed data
./GMX_fastqc.sh $OUT/qc_trimmed $OUT/reads --threads 8

# 4. Call variants  (~20 h — run it overnight)
./GMX_variant_calling.sh $OUT/reads/MYSAMPLE_R1_paired.fastq.gz \
                         $OUT/reads/MYSAMPLE_R2_paired.fastq.gz \
                         MYSAMPLE $OUT $REF 16

# 5. Annotate  (~3 h)
./GMX_annotation.sh $OUT/results $REF $FUNC

# 6. Clinical report  (~5 min)
./GMX_clinical_variants.sh $OUT/results '' $OUT/aligned_reads MYSAMPLE
```

Every script reports its version with `--version`, and its usage when run with no arguments. Each one checks its inputs first and stops
immediately with a plain message if something is missing — you will not lose hours to a typo.

### Your FASTQ files

| | |
|---|---|
| Extensions | `.fastq.gz` `.fq.gz` `.fastq.bz2` `.fq.bz2` `.fastq` `.fq`, any case |
| Pair naming | `_R1/_R2`, `_R1_001/_R2_001`, `_1/_2`, `.R1./.R2.` |
| Layout | point step 1 at a folder — subfolders are searched, so per-lane deliveries work |
| Platform | detected and reported (Illumina, MGI/DNBSEQ) |

Step 1 checks every file before anything else runs: record structure, sequence and quality
lengths, quality encoding (phred33 vs phred64), and — with `--verify` — full decompression
integrity. A truncated or non-FASTQ file stops the run then and there.

It then tells you what to pass to step 2: encoding, adapter contamination, and the MINLEN that
suits your read length.

---

## Reading your report

Variants are sorted into tiers. Read them top-down; stop when you have what you need.

| Tier | Meaning | What to do |
|---|---|---|
| **4 Secondary findings** | pathogenic variant in an ACMG SF v3.3 gene, meeting that gene's reporting rule | take to a clinician |
| **1 Known pathogenic** | ClinVar pathogenic / likely pathogenic | check the frequency column first — a "pathogenic" variant carried by 2% of people usually is not |
| **5 Carrier** | one pathogenic copy in a recessive gene | matters for family planning, usually not for your own health |
| **2 Rare loss of function** | gene switched off, rare in the population | candidates, not findings; check the quality flags |
| **3 Predicted damaging** | rare, and the predictors agree it is harmful | candidates |
| **6 For review** | conflicting ClinVar entries, rule failures, variants that turned out to be common | nothing is hidden — it is all here |

Each row carries: why it is listed, evidence codes, ACMG/AMP class and points, ClinVar
classification **and** review stars, gnomAD v4.1 frequency, inheritance, gene-disease validity,
zygosity, and quality flags.

**Two columns decide most questions.**
*ClinVar stars*: 0 stars is one anonymous submission, 4 stars is a practice guideline.
*gnomAD frequency*: if a variant is common, it does not cause a rare disease, whatever its label says.

### Sample acceptance

Before anything is classified, the sample must meet minimum quality:

| Check | Threshold |
|---|---|
| Mean genome coverage | ≥ 20x |
| Genome at ≥ 20x | ≥ 80% |
| Contamination (FREEMIX) | ≤ 3% |

Below these, the report carries a failure banner and **no variant is classified**. This is
deliberate: at low coverage the only variants that pass the depth filter are pile-ups in
collapsed repeats, which look like convincing loss-of-function findings but are artefacts.
Variants are still listed so you can see what was found — they are simply not interpreted.

Individual variants whose depth far exceeds the genome mean are flagged `DP_outlier`, for the
same reason, even in a sample that passes.

### How variants are classified

Automated ACMG/AMP scoring, ClinGen point system: pathogenic ≥ 10, likely pathogenic 6–9,
VUS 0–5, likely benign −1…−6, benign ≤ −7.

Applied: PVS1 (aware of nonsense-mediated decay escape), PM2_Supporting, PM4, PP2, PP3/BP4
(calibrated REVEL and SpliceAI), BA1, BS1, BP7.

Not applied, because no automated pipeline has the evidence: PS1–PS4, PM1, PM3, PM5, PM6, PP1,
BS2–BS4, BP1–BP3, BP5. ClinVar is shown next to the classification but never scored into it.

A variant failing quality review is not classified at all — it is marked
`Not_classified(low_quality_call)` with its flags, to be confirmed or dismissed first.

**A VUS here means not enough automated evidence. It does not mean the variant is harmless.**
The report lists every disagreement between this classification and ClinVar, so you can see
exactly where automation stops.

---

## Runtime and disk, per stage

Measured on 24 cores, 62 GB RAM, NVMe, one 30x human genome.

| Step | Time | Disk produced |
|---|---|---|
| 1, 3 Read QC | 10–30 min | < 1 GB |
| 2 Trim | 2–4 h | ~100 GB |
| 4 Call variants | ~20 h | ~110 GB |
| 5 Annotate | ~3 h | ~10 GB |
| 6 Clinical report | 30 s – 10 min | < 1 GB |

Step 6 is slow the first time only, while it looks up gnomAD v4.1 frequencies; results are
cached, so later runs take seconds. Use `--offline` to skip the lookup entirely.

---

## Troubleshooting

**`no FASTQ files found`** — check the extension is recognised, and that you pointed at the
folder, not a parent several levels up.

**Step 1 reports FAIL and stops** — the file is truncated or is not FASTQ. Re-download it.
Confirm with `--verify`.

**`gnomAD lookup failed for chrN`** — a dropped connection. Just rerun step 6; the cache keeps
what already succeeded and only the missing variants are retried.

**Funcotator fails or hangs with no network** — gnomAD 2.1 is read from cloud storage during
annotation. Step 5 needs internet.

**Out of disk during step 4** — variant calling needs ~110 GB free beyond the trimmed reads.

**GATK out of memory** — lower `--threads`, or give the JVM more heap via `_JAVA_OPTIONS`.

**Everything is `Not_classified(low_quality_call)`** — read the QC section of the report first.
Low coverage or contamination will do this, and correctly so.

---

## What this does not do

Not for exome or targeted panels · no copy-number or structural variants · no repeat expansions · no methylation · no
pharmacogenomics · mitochondrial variants are called with a diploid model and heteroplasmy is
not measured · no phenotype input, so candidates are not ranked against your symptoms ·
SpliceAI covers single-nucleotide variants only · BS1 uses global thresholds, not
gene-specific ones · short reads cannot resolve homologous or repetitive regions.

Several of these need a different technology, not just more code.

---

## Roadmap

v0.1.0 covers germline short variants end to end. These are planned, not missing pieces of it —
each needs its own tools and its own validation.

| Planned | What it adds | Why it is separate |
|---|---|---|
| Ensembl VEP annotation | a second annotation engine alongside Funcotator | different engine; needs a side-by-side comparison on the same sample before it can replace anything |
| Copy-number and structural variants | large deletions, duplications, inversions, translocations | different callers (Manta, CNVnator/Delly) and a different evidence model |
| Mitochondrial module | heteroplasmy, mtDNA variants called properly | needs Mutect2 in mitochondrial mode; the diploid caller used here cannot measure heteroplasmy |
| Pharmacogenomics | drug-response genotypes and star alleles | separate knowledge base (PharmGKB/CPIC) and haplotype calling |
| Phenotype-driven ranking | candidates ranked against the patient's HPO terms | needs phenotype input and a ranking model |
| ACMG carrier panel | the 113-gene carrier screening list | gene-list resource, then a dedicated carrier section |
| PS1 / PM5 / PM3 evidence | more ACMG codes reaching a classification | PS1/PM5 need a ClinVar protein-change index; PM3 needs variant phasing |
| Repeat expansions | STR disorders (Huntington, fragile X, ataxias) | dedicated callers (ExpansionHunter); short reads alone are limited |

Exome and targeted panels are out of scope by design — see the note at the top.

## Data sources

GRCh38 with MANE Select transcripts · Gencode v43 · ClinVar · gnomAD v2.1 and v4.1 (grpmax
FAF95) · REVEL · AlphaMissense · SpliceAI · ACMG SF v3.3 · GenCC (CC0 1.0) · gnomAD gene
constraint.

You are responsible for the licence of every source you download. **REVEL and SpliceAI are not
licensed for commercial use.**

## Licence

Proprietary. © 2026 GeneMolX AI LLC. All rights reserved. See [LICENSE](LICENSE).

Free for personal, non-commercial use. No redistribution, no commercial use, no medical use.
