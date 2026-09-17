#!/bin/bash
# © 2026 GeneMolX AI LLC. All rights reserved.
# 1500 N GRANT ST STE N DENVER, CO, 80203, USA

set -euo pipefail

# Activate genemolx-wgs environment (required for GATK, BWA, Samtools, FastQC)
GMX_VERSION="0.1.0"
case "${1:-}" in
	-V|--version) echo "GMX_variant_calling.sh ${GMX_VERSION} - GeneMolX WGS pipeline"; exit 0 ;;
esac

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate genemolx-wgs

# Script to call germline variants in a human WGS paired-end reads 
# Following GATK4 best practices workflow
# Implements the GeneMolX methodology for any human whole-genome sample.
# NOTE: all inputs are passed as arguments - review your read file names and
# paths before running. Nothing is hard-coded to a particular sample.

###################################################### VARIANT CALLING STEPS ####################################################################

# Usage: GMX_variant_calling.sh <R1_paired.fastq.gz> <R2_paired.fastq.gz> <sample_id> <output_dir> <reference_dir> [threads] [platform]
if [ $# -lt 5 ]; then
	echo "Usage: GMX_variant_calling.sh <R1_paired.fq.gz> <R2_paired.fq.gz> <sample_id> <output_dir> <reference_dir> [threads] [platform]"
	echo "  reference_dir must contain hg38.fa and the BQSR known-sites files"
	echo "  (build it with: GMX_download_references.sh <root> --reference --known-sites)"
	exit 1
fi

R1="$1"
R2="$2"
sample="$3"
output_dir="$4"
reference_dir="$5"
threads="${6:-16}"
platform="${7:-ILLUMINA}"

# Directories
ref="${reference_dir}/hg38.fa"
known_sites="${reference_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
known_indels="${reference_dir}/Homo_sapiens_assembly38.known_indels.vcf.gz"
known_mills="${reference_dir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
aligned_reads="${output_dir}/aligned_reads"
results="${output_dir}/results"
data="${output_dir}/data"
tmp="${output_dir}/tmp"
verifybamid2_svd="$(ls -d "${CONDA_PREFIX}"/share/verifybamid2-*/resource | head -1)/1000g.phase3.100k.b38.vcf.gz.dat"

mkdir -p "${aligned_reads}" "${results}" "${data}"


# Verify inputs before any long-running step starts
for f in "${R1}" "${R2}"; do
	[ -f "${f}" ] || { echo "ERROR: input FASTQ not found: ${f}"; exit 1; }
done
[ -f "${ref}" ] || { echo "ERROR: reference not found: ${ref} - run GMX_download_references.sh --reference"; exit 1; }
[ -f "${ref}.fai" ] && [ -f "${ref%.fa}.dict" ] || echo "NOTE: reference index or dictionary missing; GATK may fail"
for k in "${known_sites}" "${known_indels}" "${known_mills}"; do
	[ -f "${k}" ] || { echo "ERROR: BQSR known-sites file not found: ${k} - run GMX_download_references.sh --known-sites"; exit 1; }
done

# -------------------
# STEP 1: QC - Run fastqc --nogroup 
# -------------------

echo "STEP 1: QC - Run fastqc --nogroup"

# Disable assistive technologies and force Java into headless mode
export _JAVA_OPTIONS='-Djava.awt.headless=true -Djavax.accessibility.assistive_technologies=" "'


# QC for R1 and R2
# Read QC is a separate step - see GMX_fastqc.sh

# No trimming required, quality looks okay.



# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------

echo "STEP 2: Map to reference using BWA-MEM"


# BWA index reference 
#bwa index ${ref}

# BWA alignment
# NOTE: If you want to modify read group (RG), you can do it here manually if needed.
# ID, SM, and PL fields must match sample name conventions
bwa mem -K 100000000 -Y -t "${threads}" -R "@RG\tID:${sample}\tLB:${sample}\tPL:${platform}\tSM:${sample}" "${ref}" "${R1}" "${R2}" | samtools view -@ 4 -b -o "${aligned_reads}/sample_aligned.bam" -



# -----------------------------------------
# STEP 3: Mark Duplicates and Sort - GATK4
# -----------------------------------------

echo "STEP 3: Mark Duplicates and Sort - GATK4"


mkdir -p "${tmp}"
gatk MarkDuplicatesSpark -I "${aligned_reads}/sample_aligned.bam" -O "${aligned_reads}/sample_sorted_dedup_reads.bam" -M "${aligned_reads}/dedup_metrics.txt" --tmp-dir "${tmp}"
if [ $? -ne 0 ]; then echo "ERROR: MarkDuplicatesSpark failed"; exit 1; fi
rm "${aligned_reads}/sample_aligned.bam"



# ----------------------------------
# STEP 4: Base quality recalibration
# ----------------------------------

echo "STEP 4: Base quality recalibration"


# 1. build the model
gatk BaseRecalibrator -I "${aligned_reads}/sample_sorted_dedup_reads.bam" -R "${ref}" --known-sites "${known_sites}" --known-sites "${known_indels}" --known-sites "${known_mills}" -O "${data}/recal_data.table"

# 2. Apply the model to adjust the base quality scores
gatk ApplyBQSR -I "${aligned_reads}/sample_sorted_dedup_reads.bam" -R "${ref}" --bqsr-recal-file "${data}/recal_data.table" -O "${aligned_reads}/sample_sorted_dedup_bqsr_reads.bam"

# 3. Coverage depth QC
gatk CollectWgsMetrics \
    -I "${aligned_reads}/sample_sorted_dedup_bqsr_reads.bam" \
    -R "${ref}" \
    -O "${aligned_reads}/wgs_metrics.txt"

# 4. Contamination check
verifybamid2 \
    --SVDPrefix "${verifybamid2_svd}" \
    --Reference "${ref}" \
    --BamFile "${aligned_reads}/sample_sorted_dedup_bqsr_reads.bam" \
    --Output "${aligned_reads}/verifybamid2"
freemix="$(awk 'NR==2{print $7}' "${aligned_reads}/verifybamid2.selfSM")"
echo "FREEMIX (contamination estimate): ${freemix}"



# -----------------------------------------------
# STEP 5: Collect Alignment & Insert Size Metrics
# -----------------------------------------------

echo "STEP 5: Collect Alignment & Insert Size Metrics"


gatk CollectAlignmentSummaryMetrics R="${ref}" I="${aligned_reads}/sample_sorted_dedup_bqsr_reads.bam" O="${aligned_reads}/alignment_metrics.txt"
gatk CollectInsertSizeMetrics INPUT="${aligned_reads}/sample_sorted_dedup_bqsr_reads.bam" OUTPUT="${aligned_reads}/insert_size_metrics.txt" HISTOGRAM_FILE="${aligned_reads}/insert_size_histogram.pdf"



# ----------------------------------------------
# STEP 6: Call Variants - GATK HaplotypeCaller
# ----------------------------------------------

echo "STEP 6: Call Variants - GATK HaplotypeCaller"

gatk HaplotypeCaller -R "${ref}" -I "${aligned_reads}/sample_sorted_dedup_bqsr_reads.bam" -O "${results}/raw_variants.g.vcf.gz" -ERC GVCF --contamination-fraction-to-filter "${freemix}" --native-pair-hmm-threads 8

# Genotype the GVCF

gatk GenotypeGVCFs -R "${ref}" -V "${results}/raw_variants.g.vcf.gz" -O "${results}/raw_variants.vcf"

# Extract SNPs & INDELS

gatk SelectVariants -R "${ref}" -V "${results}/raw_variants.vcf" --select-type SNP -O "${results}/raw_snps.vcf"
gatk SelectVariants -R "${ref}" -V "${results}/raw_variants.vcf" --select-type INDEL -O "${results}/raw_indels.vcf"
