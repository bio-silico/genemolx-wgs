#!/bin/bash
# © 2026 GeneMolX AI LLC. All rights reserved.
# 1500 N GRANT ST STE N DENVER, CO, 80203, USA

set -euo pipefail

# Filters and annotates germline variants from any human whole-genome sample
# Following GATK4 best practices workflow
# Input: raw_snps.vcf and raw_indels.vcf from GMX_variant_calling.sh
# Output: curated_snps_clinical.tsv and curated_indels_clinical.tsv (51 columns)

GMX_VERSION="0.2.0"
case "${1:-}" in
	-V|--version) echo "GMX_annotation.sh ${GMX_VERSION} - GeneMolX WGS pipeline"; exit 0 ;;
esac

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate genemolx-wgs

echo "Started: $(date)"

# Usage: GMX_annotation.sh <results_dir> <reference_dir> <funcotator_data_sources_dir> [acmg_gene_list]
if [ $# -lt 3 ]; then
	echo "Usage: GMX_annotation.sh <results_dir> <reference_dir> <funcotator_data_sources_dir> [acmg_gene_list]"
	echo "  results_dir must contain raw_snps.vcf and raw_indels.vcf from GMX_variant_calling.sh"
	exit 1
fi

results="$1"
reference_dir="$2"
funcotator_data="$3"
acmg_list="${4:-$(dirname "$0")/resources/acmg/ACMG_SF_v3.3_gene_list.tsv}"

# Directories
ref="${reference_dir}/hg38.fa"

# Verify inputs before any long-running step starts
[ -d "${results}" ] || { echo "ERROR: results directory not found: ${results}"; exit 1; }
for v in raw_snps.vcf raw_indels.vcf; do
	[ -f "${results}/${v}" ] || { echo "ERROR: ${v} not found in ${results} - run GMX_variant_calling.sh first"; exit 1; }
done
[ -f "${ref}" ] || { echo "ERROR: reference not found: ${ref}"; exit 1; }
[ -d "${funcotator_data}" ] || { echo "ERROR: Funcotator data sources not found: ${funcotator_data}"; exit 1; }
[ -f "${acmg_list}" ] || { echo "ERROR: ACMG gene list not found: ${acmg_list}"; exit 1; }

# -------------------
# STEP 1: Filter SNPs
# -------------------

echo "STEP 1: Filter SNPs"

gatk VariantFiltration \
	-R "${ref}" \
	-V "${results}/raw_snps.vcf" \
	-O "${results}/filtered_snps.vcf" \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 3.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter" \
	--set-filtered-genotype-to-no-call

if [ $? -ne 0 ]; then echo "ERROR: SNP VariantFiltration failed"; exit 1; fi

# -------------------
# STEP 2: Filter INDELs
# -------------------

echo "STEP 2: Filter INDELs"

gatk VariantFiltration \
	-R "${ref}" \
	-V "${results}/raw_indels.vcf" \
	-O "${results}/filtered_indels.vcf" \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -20.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter" \
	--set-filtered-genotype-to-no-call

if [ $? -ne 0 ]; then echo "ERROR: INDEL VariantFiltration failed"; exit 1; fi

# -------------------
# STEP 3: Select PASS variants
# -------------------

echo "STEP 3: Select PASS variants"

gatk SelectVariants \
	--exclude-filtered \
	-V "${results}/filtered_snps.vcf" \
	-O "${results}/analysis-ready-snps.vcf"

gatk SelectVariants \
	--exclude-filtered \
	-V "${results}/filtered_indels.vcf" \
	-O "${results}/analysis-ready-indels.vcf"

if [ $? -ne 0 ]; then echo "ERROR: SelectVariants failed"; exit 1; fi

# -------------------
# STEP 4: Exclude genotype-level filter failures
# -------------------

echo "STEP 4: Exclude genotype-level filter failures"

gatk SelectVariants \
	--exclude-non-variants \
	-V "${results}/analysis-ready-snps.vcf" \
	-O "${results}/analysis-ready-snps-filteredGT.vcf"

gatk SelectVariants \
	--exclude-non-variants \
	-V "${results}/analysis-ready-indels.vcf" \
	-O "${results}/analysis-ready-indels-filteredGT.vcf"

# -------------------
# STEP 5: Annotate with Funcotator
# -------------------

echo "STEP 5: Normalize variants (split multi-allelic sites, left-align, remove * alleles)"

for type in snps indels; do
	bcftools norm -m -both -f "${ref}" "${results}/analysis-ready-${type}-filteredGT.vcf" | \
		bcftools view -e 'ALT="*"' | \
		bcftools sort -T "${results}/bcftools_sort." -Oz -o "${results}/analysis-ready-${type}-filteredGT-norm.vcf.gz"
	tabix -p vcf "${results}/analysis-ready-${type}-filteredGT-norm.vcf.gz"
done

echo "STEP 5: Annotate SNPs with Funcotator"

gatk Funcotator \
	--variant "${results}/analysis-ready-snps-filteredGT-norm.vcf.gz" \
	--reference "${ref}" \
	--ref-version hg38 \
	--data-sources-path "${funcotator_data}" \
	--prefer-mane-transcripts \
	--output "${results}/analysis-ready-snps-filteredGT-functotated.vcf" \
	--output-file-format VCF

if [ $? -ne 0 ]; then echo "ERROR: SNP Funcotator failed"; exit 1; fi

echo "STEP 5: Annotate INDELs with Funcotator"

gatk Funcotator \
	--variant "${results}/analysis-ready-indels-filteredGT-norm.vcf.gz" \
	--reference "${ref}" \
	--ref-version hg38 \
	--data-sources-path "${funcotator_data}" \
	--prefer-mane-transcripts \
	--output "${results}/analysis-ready-indels-filteredGT-functotated.vcf" \
	--output-file-format VCF

if [ $? -ne 0 ]; then echo "ERROR: INDEL Funcotator failed"; exit 1; fi

# -------------------
# STEP 6: Convert to tables
# -------------------

echo "STEP 6: Convert to tables"

gatk VariantsToTable \
	-V "${results}/analysis-ready-snps-filteredGT-functotated.vcf" \
	-F AC -F AN -F DP -F AF -F FUNCOTATION \
	-O "${results}/output_snps.table"

gatk VariantsToTable \
	-V "${results}/analysis-ready-indels-filteredGT-functotated.vcf" \
	-F AC -F AN -F DP -F AF -F FUNCOTATION \
	-O "${results}/output_indels.table"

if [ $? -ne 0 ]; then echo "ERROR: VariantsToTable failed"; exit 1; fi

# -------------------
# STEP 7: Extract curated clinical variant tables
# -------------------
# Extracts: VCF fields (coordinates, quality, genotype) + selected Funcotation fields
# Output: Human/AI-readable TSV with clinically relevant columns
# Funcotation fields are selected by NAME from the VCF header (##INFO=<ID=FUNCOTATION ...>)

echo "STEP 7: Extract curated clinical variant tables"

# Create reusable awk script for Funcotation extraction (mawk-compatible)
AWK_SCRIPT=$(mktemp /tmp/funcotation_extract.XXXXXX.awk)
cat > "${AWK_SCRIPT}" << 'AWKEOF'
function val(name) { return ((name in idx) && (idx[name] in f)) ? f[idx[name]] : "." }
function ival(name) { return (name in iv) ? iv[name] : "." }
BEGIN {
    OFS="\t"
    # ACMG SF gene list from resources file (Gene = column 4)
    while ((getline line < acmg_file) > 0) {
        if (++nl == 1) continue
        split(line, a, "\t"); acmg[a[4]] = 1
    }
    close(acmg_file)
    # Existing columns (same order as before)
    nsel = split("Gencode_43_hugoSymbol Gencode_43_variantClassification Gencode_43_variantType Gencode_43_genomeChange Gencode_43_annotationTranscript Gencode_43_transcriptExon Gencode_43_cDnaChange Gencode_43_codonChange Gencode_43_proteinChange ACMGLMMLof_LOF_Mechanism ACMGLMMLof_Mode_of_Inheritance ACMG_recommendation_Disease_Name ClinVar_VCF_CLNSIG ClinVar_VCF_CLNREVSTAT ClinVar_VCF_CLNDN ClinVar_VCF_CLNHGVS ClinVar_VCF_RS LMMKnown_LMM_FLAGGED gnomAD_exome_AF gnomAD_exome_AF_popmax gnomAD_genome_AF gnomAD_genome_AF_popmax Gencode_43_otherTranscripts", sel, " ")
    # New columns (appended at the end)
    nnew = split("gnomAD_v4_1_AF_joint gnomAD_v4_1_fafmax_faf95_max_joint AlphaMissense_am_pathogenicity AlphaMissense_am_class REVEL_REVEL SpliceAI_DS_max", newf, " ")
}
/^##INFO=<ID=FUNCOTATION/ {
    h = $0; sub(/.*Funcotation fields are: /, "", h); sub(/">$/, "", h)
    n = split(h, names, "|"); for (i = 1; i <= n; i++) idx[names[i]] = i
    next
}
/^#/ { next }
{
    chrom=$1; pos=$2; id=$3; ref=$4; alt=$5; qual=$6
    nf=split($9, fmt, ":")
    split($10, samp, ":")
    gt=""; dp=""; ad=""; gq=""
    pgt=""; pid=""; ps=""
    for (i=1; i<=nf; i++) {
        if (fmt[i]=="GT") gt=samp[i]
        if (fmt[i]=="DP") dp=samp[i]
        if (fmt[i]=="AD") ad=samp[i]
        if (fmt[i]=="GQ") gq=samp[i]
        # physical phasing, emitted by HaplotypeCaller only where it could phase
        if (fmt[i]=="PGT") pgt=samp[i]
        if (fmt[i]=="PID") pid=samp[i]
        if (fmt[i]=="PS")  ps=samp[i]
    }
    # VAF: sites are bi-allelic after normalization (AD = REF,ALT)
    vaf="."
    if (ad != "" && ad != ".") {
        split(ad, adv, ",")
        total=adv[1]+adv[2]
        if (total>0) vaf=sprintf("%.3f", adv[2]/total)
    }
    fstr=$0
    sub(/.*FUNCOTATION=\[/, "", fstr)
    sub(/].*/, "", fstr)
    split(fstr, f, "|")
    # GATK quality annotations from INFO (FUNCOTATION block removed first)
    info=$8
    sub(/FUNCOTATION=\[[^]]*\]/, "", info)
    split("", iv)
    ni=split(info, ia, ";")
    for (i=1; i<=ni; i++) { p=index(ia[i], "="); if (p>0) iv[substr(ia[i],1,p-1)]=substr(ia[i],p+1) }
    out = chrom OFS pos OFS id OFS ref OFS alt OFS qual OFS gt OFS dp OFS ad OFS gq OFS vaf
    for (i = 1; i <= nsel; i++) out = out OFS val(sel[i])
    # ACMG_SF: True if gene is in ACMG SF list (resources file)
    acmg_sf = (val("Gencode_43_hugoSymbol") in acmg) ? "True" : "False"
    # Splice_Region: intronic offset within canonical (+/-1,2) or extended splice region
    # (donor +3..+8, acceptor -3..-12); requires a coding position before the offset (excludes 5'UTR c.-N)
    cdna = val("Gencode_43_cDnaChange")
    splice_region = "False"
    if (cdna ~ /[0-9]\+[1-8][^0-9]/ || cdna ~ /[0-9]-[1-9][^0-9]/ || cdna ~ /[0-9]-1[0-2][^0-9]/)
        splice_region = "True"
    out = out OFS acmg_sf OFS splice_region
    for (i = 1; i <= nnew; i++) out = out OFS val(newf[i])
    out = out OFS ival("QD") OFS ival("FS") OFS ival("SOR") OFS ival("MQ") OFS ival("MQRankSum") OFS ival("ReadPosRankSum")
    # phase set: needed to tell in trans from in cis for ACMG PM3
    out = out OFS (pgt=="" ? "." : pgt) OFS (pid=="" ? "." : pid) OFS (ps=="" ? "." : ps)
    print out
}
AWKEOF

HEADER="CHROM\tPOS\tID\tREF\tALT\tQUAL\tGT\tDP\tAD\tGQ\tVAF\tGene\tVariant_Classification\tVariant_Type\tGenome_Change\tTranscript\tExon\tcDNA_Change\tCodon_Change\tProtein_Change\tLOF_Mechanism\tMode_of_Inheritance\tACMG_Disease\tClinVar_Significance\tClinVar_Review_Status\tClinVar_Disease\tClinVar_HGVS\tdbSNP_RS\tLMM_Flagged\tgnomAD_Exome_AF\tgnomAD_Exome_AF_PopMax\tgnomAD_Genome_AF\tgnomAD_Genome_AF_PopMax\tOther_Transcripts\tACMG_SF\tSplice_Region\tgnomAD_v4.1_Joint_AF\tgnomAD_v4.1_Grpmax_FAF95\tAlphaMissense_Score\tAlphaMissense_Class\tREVEL_Score\tSpliceAI_DS_Max\tQD\tFS\tSOR\tMQ\tMQRankSum\tReadPosRankSum\tPGT\tPID\tPS"

# SNPs
echo "  Extracting curated SNP table..."
printf "${HEADER}\n" > "${results}/curated_snps_clinical.tsv"
mawk -F'\t' -v acmg_file="${acmg_list}" -f "${AWK_SCRIPT}" "${results}/analysis-ready-snps-filteredGT-functotated.vcf" >> "${results}/curated_snps_clinical.tsv"
snp_count=$(( $(wc -l < "${results}/curated_snps_clinical.tsv") - 1 ))
echo "  SNPs: ${snp_count} variants"

# INDELs
echo "  Extracting curated INDEL table..."
printf "${HEADER}\n" > "${results}/curated_indels_clinical.tsv"
mawk -F'\t' -v acmg_file="${acmg_list}" -f "${AWK_SCRIPT}" "${results}/analysis-ready-indels-filteredGT-functotated.vcf" >> "${results}/curated_indels_clinical.tsv"
indel_count=$(( $(wc -l < "${results}/curated_indels_clinical.tsv") - 1 ))
echo "  INDELs: ${indel_count} variants"

rm -f "${AWK_SCRIPT}"

echo "Variant filtering, annotation, and table extraction completed!"
echo "Finished: $(date)"
