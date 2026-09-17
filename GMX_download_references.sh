#!/bin/bash
# © 2026 GeneMolX AI LLC. All rights reserved.
# 1500 N GRANT ST STE N DENVER, CO, 80203, USA

set -euo pipefail

# Script to download and prepare reference data for the GeneMolX WGS pipeline
# Output layout under <reference_root>:
#   hg38/                                          -> reference_dir for GMX_variant_calling.sh
#   funcotator_dataSources.v1.8.hg38.20230908g/    -> funcotator_data_sources_dir for GMX_annotation.sh
#
# Usage: GMX_download_references.sh <reference_root> [--all] [--reference] [--known-sites] [--funcotator]
#                                   [--clinvar] [--alphamissense] [--revel] [--spliceai] [--gencc] [--constraint] [--exon-structure]
#   --all runs every section
#   --clinvar, --alphamissense, --revel, --spliceai need the Funcotator bundle (--funcotator) first
#   --gencc writes a gene-level inheritance / gene-disease validity table into resources/gencc/
#     (GenCC submissions export, CC0 1.0; used by GMX_clinical_variants.sh for tiers 4 and 5)
#
# SpliceAI: public Ensembl precomputed masked scores for all SNVs in MANE Select v1.4 transcripts (GRCh38).
#   Ensembl does not provide SpliceAI scores for indels (indels only from Illumina BaseSpace, login required).

GMX_VERSION="0.1.0"
case "${1:-}" in
	-V|--version) echo "GMX_download_references.sh ${GMX_VERSION} - GeneMolX WGS pipeline"; exit 0 ;;
esac

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate genemolx-wgs

echo "Started: $(date)"

root="${1:?Usage: GMX_download_references.sh <reference_root> [sections]}"
shift
do_ref=0; do_known=0; do_func=0; do_clinvar=0; do_am=0; do_revel=0; do_spliceai=0; do_gencc=0; do_constraint=0; do_exons=0
[ $# -eq 0 ] && { echo "ERROR: no section selected"; exit 1; }
while [ $# -gt 0 ]; do
	case "$1" in
		--all) do_ref=1; do_known=1; do_func=1; do_clinvar=1; do_am=1; do_revel=1; do_spliceai=1; do_gencc=1; do_constraint=1; do_exons=1 ;;
		--reference) do_ref=1 ;;
		--known-sites) do_known=1 ;;
		--funcotator) do_func=1 ;;
		--clinvar) do_clinvar=1 ;;
		--alphamissense) do_am=1 ;;
		--revel) do_revel=1 ;;
		--spliceai) do_spliceai=1 ;;
		--gencc) do_gencc=1 ;;
		--constraint) do_constraint=1 ;;
		--exon-structure) do_exons=1 ;;
		*) echo "ERROR: unknown option $1"; exit 1 ;;
	esac
	shift
done

# Directories
ref_dir="${root}/hg38"
bundle="${root}/funcotator_dataSources.v1.8.hg38.20230908g"
dl="${root}/downloads"
tmp="${root}/tmp"
mkdir -p "${ref_dir}" "${dl}" "${tmp}"

# Helpers
fetch() { curl -fL --retry 5 --retry-delay 10 -C - -o "$2" "$1"; }
check_md5() { echo "$2  $1" | md5sum -c -; }
check_sha256() { echo "$2  $1" | sha256sum -c -; }
contig_header() { [ -s "${ref_dir}/hg38.fa.fai" ] && awk '{print "##contig=<ID=" $1 ",length=" $2 ">"}' "${ref_dir}/hg38.fa.fai" || true; }
write_vcf_config() {  # <config_path> <name> <version> <src_file> <origin_location>
	cat > "$1" << EOF
name = $2
version = $3
src_file = $4
origin_location = $5
preprocessing_script = GMX_download_references.sh

# Whether this data source is for the b37 reference.
# Required and defaults to false.
isB37DataSource = false

# Supported types:
# simpleXSV    -- Arbitrary separated value table (e.g. CSV), keyed off Gene Name OR Transcript ID
# locatableXSV -- Arbitrary separated value table (e.g. CSV), keyed off a genome location
# gencode      -- Custom datasource class for GENCODE
#	cosmic       -- Custom datasource class for COSMIC
# vcf          -- Custom datasource class for Variant Call Format (VCF) files
type = vcf

# Required field for GENCODE files.
# Path to the FASTA file from which to load the sequences for GENCODE transcripts:
gencode_fasta_path =

# Required field for GENCODE files.
# NCBI build version (either hg19 or hg38):
ncbi_build_version =

# Required field for simpleXSV files.
# Valid values:
#     GENE_NAME
#     TRANSCRIPT_ID
xsv_key = GENE_NAME

# Required field for simpleXSV files.
# The 0-based index of the column containing the key on which to match
xsv_key_column = 0

# Required field for simpleXSV AND locatableXSV files.
# The delimiter by which to split the XSV file into columns.
xsv_delimiter = ,

# Required field for simpleXSV files.
# Whether to permissively match the number of columns in the header and data rows
# Valid values:
#     true
#     false
xsv_permissive_cols = true

# Required field for locatableXSV files.
# The 0-based index of the column containing the contig for each row
contig_column =

# Required field for locatableXSV files.
# The 0-based index of the column containing the start position for each row
start_column =

# Required field for locatableXSV files.
# The 0-based index of the column containing the end position for each row
end_column =
EOF
}

# -------------------
# SECTION 1: Reference genome (UCSC hg38) + indexes
# -------------------

if [ "${do_ref}" -eq 1 ]; then
	echo "SECTION 1: Reference genome hg38"
	base="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips"
	fetch "${base}/hg38.fa.gz" "${dl}/hg38.fa.gz"
	fetch "${base}/md5sum.txt" "${dl}/hg38_md5sum.txt"
	check_md5 "${dl}/hg38.fa.gz" "$(awk '$2 == "hg38.fa.gz" {print $1}' "${dl}/hg38_md5sum.txt")"
	zcat "${dl}/hg38.fa.gz" > "${ref_dir}/hg38.fa"
	# Validated GeneMolX reference (hg38.fa uncompressed)
	check_md5 "${ref_dir}/hg38.fa" "b2aee9f885accc00531e59c4736bee63" || echo "WARNING: hg38.fa differs from the validated GeneMolX reference"
	samtools faidx "${ref_dir}/hg38.fa"
	gatk CreateSequenceDictionary -R "${ref_dir}/hg38.fa" -O "${ref_dir}/hg38.dict"
	bwa index "${ref_dir}/hg38.fa"
fi

# -------------------
# SECTION 2: Known sites for BQSR (GATK resource bundle hg38 v0)
# -------------------

if [ "${do_known}" -eq 1 ]; then
	echo "SECTION 2: Known sites (dbSNP138, Mills and 1000G, known indels)"
	base="https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"
	fetch "${base}/Homo_sapiens_assembly38.dbsnp138.vcf" "${dl}/Homo_sapiens_assembly38.dbsnp138.vcf"
	check_md5 "${dl}/Homo_sapiens_assembly38.dbsnp138.vcf" "f7e1ef5c1830bfb33675b9c7cbaa4868"
	bgzip -c "${dl}/Homo_sapiens_assembly38.dbsnp138.vcf" > "${ref_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
	tabix -p vcf "${ref_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
	[ "$(bcftools index -n "${ref_dir}/Homo_sapiens_assembly38.dbsnp138.vcf.gz")" -eq 60691395 ] || { echo "ERROR: dbSNP record count"; exit 1; }

	for f in Mills_and_1000G_gold_standard.indels.hg38.vcf.gz:2e02696032dcfe95ff0324f4a13508e3 \
	         Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi:4c807e2cbe0752c0c44ac82ff3b52025 \
	         Homo_sapiens_assembly38.known_indels.vcf.gz:14cc588a271951ac1806f9be895fb51f \
	         Homo_sapiens_assembly38.known_indels.vcf.gz.tbi:1a55fdfa6533ae5cbc70e8188e779229; do
		name="${f%%:*}"; md5="${f##*:}"
		fetch "${base}/${name}" "${ref_dir}/${name}"
		check_md5 "${ref_dir}/${name}" "${md5}"
	done
fi

# -------------------
# SECTION 3: Funcotator data sources (germline, hg38) + gnomAD 2.1 (remote)
# -------------------

if [ "${do_func}" -eq 1 ]; then
	echo "SECTION 3: Funcotator data sources v1.8"
	base="https://storage.googleapis.com/gcp-public-data--broad-references/funcotator"
	fetch "${base}/funcotator_dataSources.v1.8.hg38.20230908g.tar.gz" "${dl}/funcotator_dataSources.v1.8.hg38.20230908g.tar.gz"
	check_sha256 "${dl}/funcotator_dataSources.v1.8.hg38.20230908g.tar.gz" "a7ba51989aef8438dd8f32705152c68f9aeebed07d6b13391c5ca049f368c5cc"
	tar -xzf "${dl}/funcotator_dataSources.v1.8.hg38.20230908g.tar.gz" -C "${root}"
	[ -d "${bundle}" ] || { echo "ERROR: ${bundle} not found after extraction"; exit 1; }
	# Enable gnomAD 2.1 (exome + genome); original bucket is requester-pays -> public broad-references bucket
	tar -xzf "${bundle}/gnomAD_exome.tar.gz" -C "${bundle}"
	tar -xzf "${bundle}/gnomAD_genome.tar.gz" -C "${bundle}"
	sed -i 's#gs://broad-public-datasets/funcotator/#gs://gcp-public-data--broad-references/funcotator/#' \
		"${bundle}/gnomAD_exome/hg38/gnomAD_exome.config" "${bundle}/gnomAD_genome/hg38/gnomAD_genome.config"
fi

# -------------------
# SECTION 4: ClinVar (current weekly release, GRCh38) -> Funcotator ClinVar_VCF
# -------------------

if [ "${do_clinvar}" -eq 1 ]; then
	echo "SECTION 4: ClinVar"
	[ -d "${bundle}/clinvar/hg38" ] || { echo "ERROR: run --funcotator first"; exit 1; }
	base="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38"
	fetch "${base}/clinvar.vcf.gz" "${dl}/clinvar.vcf.gz"
	fetch "${base}/clinvar.vcf.gz.md5" "${dl}/clinvar.vcf.gz.md5"
	check_md5 "${dl}/clinvar.vcf.gz" "$(awk '{print $1}' "${dl}/clinvar.vcf.gz.md5")"
	release="$(grep -o 'clinvar_[0-9]\{8\}' "${dl}/clinvar.vcf.gz.md5" | head -1 | cut -d_ -f2)"
	out="clinvar_${release}_hg38_chr.vcf.gz"
	# 1..22, X, Y -> chr1..chrY ; MT -> chrM ; other contigs (NT_/NW_) dropped
	zcat "${dl}/clinvar.vcf.gz" | \
		awk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} $1 ~ /^([0-9]+|X|Y|MT)$/ {$1 = ($1 == "MT") ? "chrM" : "chr" $1; print}' | \
		bgzip -c > "${bundle}/clinvar/hg38/${out}"
	tabix -p vcf "${bundle}/clinvar/hg38/${out}"
	sed -i -e "s#^version = .*#version = ${release}_hg38#" -e "s#^src_file = .*#src_file = ${out}#" \
		"${bundle}/clinvar/hg38/clinvar_vcf.config"
fi

# -------------------
# SECTION 5: AlphaMissense (hg38) -> Funcotator AlphaMissense
# -------------------

if [ "${do_am}" -eq 1 ]; then
	echo "SECTION 5: AlphaMissense"
	[ -d "${bundle}" ] || { echo "ERROR: run --funcotator first"; exit 1; }
	fetch "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz" "${dl}/AlphaMissense_hg38.tsv.gz"
	mkdir -p "${bundle}/AlphaMissense/hg38"
	out="${bundle}/AlphaMissense/hg38/AlphaMissense_hg38.vcf.gz"
	{
		printf '##fileformat=VCFv4.2\n'
		printf '##INFO=<ID=am_pathogenicity,Number=1,Type=Float,Description="AlphaMissense pathogenicity score">\n'
		printf '##INFO=<ID=am_class,Number=1,Type=String,Description="AlphaMissense class">\n'
		contig_header
		printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
		zcat "${dl}/AlphaMissense_hg38.tsv.gz" | \
			awk 'BEGIN{FS=OFS="\t"} !/^#/ {print $1, $2, ".", $3, $4, ".", ".", "am_pathogenicity=" $9 ";am_class=" $10}'
	} | bcftools sort -m 2G -T "${tmp}" -Oz -o "${out}"
	tabix -p vcf "${out}"
	write_vcf_config "${bundle}/AlphaMissense/hg38/alphamissense.config" "AlphaMissense" "hg38" \
		"AlphaMissense_hg38.vcf.gz" "https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
fi

# -------------------
# SECTION 6: REVEL v1.3 -> Funcotator REVEL
# -------------------

if [ "${do_revel}" -eq 1 ]; then
	echo "SECTION 6: REVEL"
	[ -d "${bundle}" ] || { echo "ERROR: run --funcotator first"; exit 1; }
	command -v unzip > /dev/null || { echo "ERROR: unzip not found"; exit 1; }
	fetch "https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip" "${dl}/revel-v1.3_all_chromosomes.zip"
	mkdir -p "${bundle}/REVEL/hg38"
	out="${bundle}/REVEL/hg38/revel_v1.3_grch38.vcf.gz"
	{
		printf '##fileformat=VCFv4.2\n'
		printf '##INFO=<ID=REVEL,Number=1,Type=Float,Description="REVEL score v1.3 (max across transcripts)">\n'
		contig_header
		printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
		# CSV: chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid (one row per transcript)
		unzip -p "${dl}/revel-v1.3_all_chromosomes.zip" revel_with_transcript_ids | \
			awk 'BEGIN{FS=","; OFS="\t"} NR > 1 && $3 != "." {print "chr" $1, $3, $4, $5, $8}' | \
			LC_ALL=C sort -S 4G -T "${tmp}" -k1,1 -k2,2n -k3,3 -k4,4 | \
			awk 'BEGIN{FS=OFS="\t"}
			     { key = $1 FS $2 FS $3 FS $4
			       if (key == prev) { if ($5 > max) max = $5; next }
			       if (NR > 1) print c, p, ".", r, a, ".", ".", "REVEL=" max
			       prev = key; c = $1; p = $2; r = $3; a = $4; max = $5 }
			     END { if (NR > 0) print c, p, ".", r, a, ".", ".", "REVEL=" max }'
	} | bcftools sort -m 2G -T "${tmp}" -Oz -o "${out}"
	tabix -p vcf "${out}"
	write_vcf_config "${bundle}/REVEL/hg38/revel.config" "REVEL" "v1.3" \
		"revel_v1.3_grch38.vcf.gz" "https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip"
fi

# -------------------
# SECTION 7: SpliceAI (Ensembl precomputed masked SNV scores, MANE Select v1.4) -> Funcotator SpliceAI (DS_max)
# -------------------

if [ "${do_spliceai}" -eq 1 ]; then
	echo "SECTION 7: SpliceAI (Ensembl, masked SNV scores, MANE Select v1.4)"
	[ -d "${bundle}" ] || { echo "ERROR: run --funcotator first"; exit 1; }
	base="https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_plugins"
	src="spliceai_scores.masked.snv.ensembl_mane_v1.4.grch38.vcf.gz"
	fetch "${base}/${src}" "${dl}/${src}"
	# No published checksum: check BGZF integrity of the full download
	bgzip -t "${dl}/${src}"
	mkdir -p "${bundle}/SpliceAI/hg38"
	out="${bundle}/SpliceAI/hg38/spliceai_masked_snv_mane_v1.4_hg38_dsmax.vcf.gz"
	# INFO SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL (comma-separated per gene)
	# Contigs 1..22, X, Y -> chr1..chrY ; MT -> chrM
	zcat "${dl}/${src}" | awk 'BEGIN{FS=OFS="\t"}
		/^##fileformat/ {print; print "##INFO=<ID=DS_max,Number=1,Type=Float,Description=\"SpliceAI max delta score (max of DS_AG, DS_AL, DS_DG, DS_DL across genes)\">"; next}
		/^##contig=<ID=/ { sub(/ID=MT,/, "ID=chrM,"); if ($0 !~ /ID=chr/) sub(/ID=/, "ID=chr"); print; next }
		/^##/ {next}
		/^#CHROM/ {print; next}
		{ if ($1 == "MT") $1 = "chrM"; else if ($1 !~ /^chr/) $1 = "chr" $1
		  s = $8; sub(/.*SpliceAI=/, "", s); sub(/;.*/, "", s)
		  m = 0; n = split(s, e, ",")
		  for (i = 1; i <= n; i++) { split(e[i], v, "|"); for (j = 3; j <= 6; j++) if (v[j] + 0 > m) m = v[j] + 0 }
		  $8 = "DS_max=" m; print }' | bgzip -c > "${out}"
	tabix -p vcf "${out}"
	write_vcf_config "${bundle}/SpliceAI/hg38/spliceai.config" "SpliceAI" "ensembl_mane_v1.4_masked_snv" \
		"$(basename "${out}")" "${base}/${src}"
fi

# -------------------
# SECTION 8: FUNCOTATION field names (must match GMX_annotation.sh Step 7)
# -------------------

if [ -d "${bundle}" ]; then
	echo "SECTION 8: FUNCOTATION field names from VCF data sources"
	for ds in clinvar AlphaMissense REVEL SpliceAI; do
		cfg=$(ls "${bundle}/${ds}/hg38/"*.config 2> /dev/null | head -1 || true)
		[ -n "${cfg}" ] || continue
		name=$(awk -F' = ' '$1 == "name" {print $2}' "${cfg}")
		src=$(awk -F' = ' '$1 == "src_file" {print $2}' "${cfg}")
		echo "  ${name}: $(bcftools view -h "${bundle}/${ds}/hg38/${src}" | sed -n "s/^##INFO=<ID=\([^,]*\),.*/${name}_\1/p" | tr '\n' ' ')"
	done
fi

# -------------------
# SECTION 9: GenCC gene-disease validity and mode of inheritance
# -------------------
# Source: GenCC submissions export (ClinGen, Genomics England PanelApp, Orphanet,
# Ambry, Invitae, Illumina and others), CC0 1.0 Public Domain Dedication.
# Collapsed to one row per gene: strongest classification, observed modes of
# inheritance, and the disease of the strongest submission.

if [ "${do_gencc}" -eq 1 ]; then
	echo "SECTION 9: GenCC gene-disease validity"
	gencc_dir="$(dirname "$0")/resources/gencc"
	mkdir -p "${gencc_dir}"
	fetch "https://search.thegencc.org/download/action/submissions-export-tsv" "${gencc_dir}/gencc-submissions.tsv"

	mawk -F'\t' '
	function rank(c) {
		c = tolower(c)
		if (c ~ /definitive/)  return 6
		if (c ~ /strong/)      return 5
		if (c ~ /moderate/)    return 4
		if (c ~ /supportive/)  return 3
		if (c ~ /limited/)     return 2
		if (c ~ /disputed/)    return 1
		return 0        # refuted, animal model only, no known disease relationship
	}
	function moi(m) {
		m = tolower(m)
		if (m ~ /autosomal recessive/) return "AR"
		if (m ~ /autosomal dominant/)  return "AD"
		if (m ~ /x-linked recessive/)  return "XLR"
		if (m ~ /x-linked dominant/)   return "XLD"
		if (m ~ /x-linked/)            return "XL"
		if (m ~ /mitochondrial/)       return "MT"
		if (m ~ /semidominant/)        return "SD"
		if (m ~ /somatic/)             return "SOM"
		return ""
	}
	NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
	{
		g = $(h["gene_symbol"]); if (g == "") next
		gsub(/"/, "", g)
		r = rank($(h["classification_title"]))
		m = moi($(h["moi_title"]))
		d = $(h["disease_title"]); gsub(/"/, "", d)
		if (!(g in best) || r > best[g]) { best[g] = r; bestc[g] = $(h["classification_title"]); bestd[g] = d }
		# inheritance is only taken from submissions with at least limited evidence
		if (m != "" && r >= 2 && index("/" inh[g] "/", "/" m "/") == 0)
			inh[g] = (inh[g] == "" ? m : inh[g] "/" m)
	}
	END {
		OFS = "\t"
		print "Gene", "Inheritance", "Validity", "Disease"
		for (g in best) {
			c = bestc[g]; gsub(/"/, "", c)
			print g, (inh[g] == "" ? "." : inh[g]), (c == "" ? "." : c), (bestd[g] == "" ? "." : bestd[g])
		}
	}' "${gencc_dir}/gencc-submissions.tsv" | { read -r hdr; echo "${hdr}"; sort; } > "${gencc_dir}/gencc_gene_inheritance.tsv"

	echo "  genes: $(( $(wc -l < "${gencc_dir}/gencc_gene_inheritance.tsv") - 1 ))"
	echo "  ${gencc_dir}/gencc_gene_inheritance.tsv"
fi

# -------------------
# SECTION 10: gnomAD v4.1 gene constraint metrics
# -------------------
# LOEUF (lof.oe_ci.upper), missense z-score and pLI per gene, used by the
# ACMG/AMP classifier: PVS1 needs a gene where loss of function is a disease
# mechanism, PP2 needs a gene constrained against missense variation.

if [ "${do_constraint}" -eq 1 ]; then
	echo "SECTION 10: gnomAD v4.1 constraint metrics"
	con_dir="$(dirname "$0")/resources/gnomad"
	mkdir -p "${con_dir}"
	fetch "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv" \
	      "${tmp}/gnomad.v4.1.constraint_metrics.tsv"

	mawk -F'\t' '
	NR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
	{
		g = $(h["gene"]); if (g == "" || g == "NA") next
		mane = $(h["mane_select"])
		# one row per gene: prefer the MANE Select transcript
		if ((g in seen) && seen[g] == "true" && mane != "true") next
		seen[g] = mane
		loeuf[g] = $(h["lof.oe_ci.upper"]); misz[g] = $(h["mis.z_score"]); pli[g] = $(h["lof.pLI"])
	}
	END {
		OFS = "\t"
		print "Gene", "LOEUF", "Mis_Z", "pLI"
		for (g in loeuf)
			print g, (loeuf[g] == "" ? "NA" : loeuf[g]), (misz[g] == "" ? "NA" : misz[g]), (pli[g] == "" ? "NA" : pli[g])
	}' "${tmp}/gnomad.v4.1.constraint_metrics.tsv" | { read -r hdr; echo "${hdr}"; sort; } > "${con_dir}/gnomad_v4.1_constraint.tsv"

	echo "  genes: $(( $(wc -l < "${con_dir}/gnomad_v4.1_constraint.tsv") - 1 ))"
fi

# -------------------
# SECTION 11: transcript exon structure (for PVS1 / NMD escape)
# -------------------
# From the Gencode GTF already present in the Funcotator bundle. PVS1 is
# downgraded when a truncating variant escapes nonsense-mediated decay, i.e.
# it falls in the last exon or the last 50 bp of the penultimate exon.

if [ "${do_exons}" -eq 1 ]; then
	echo "SECTION 11: transcript exon structure"
	gtf=$(ls "${bundle}/gencode/hg38/"*.gtf 2> /dev/null | head -1 || true)
	if [ -z "${gtf}" ]; then
		echo "  ERROR: no Gencode GTF in ${bundle}/gencode/hg38 (run --funcotator first)"; exit 1
	fi
	ex_dir="$(dirname "$0")/resources/gencode"
	mkdir -p "${ex_dir}"

	mawk -F'\t' '
	$3 == "exon" {
		if (match($9, /transcript_id "[^"]+"/)) {
			t = substr($9, RSTART + 15, RLENGTH - 16)
			sub(/\..*$/, "", t)                       # drop the version
			n[t]++
			if (match($9, /exon_number [0-9]+/)) {
				e = substr($9, RSTART + 12, RLENGTH - 12) + 0
				if (e > maxe[t]) maxe[t] = e
				st[t "_" e] = $4; en[t "_" e] = $5
			}
			strand[t] = $7
		}
	}
	END {
		OFS = "\t"
		print "Transcript", "N_Exons", "Strand", "Penult_Start", "Penult_End", "Last_Start", "Last_End"
		for (t in maxe) {
			m = maxe[t]; p = m - 1
			print t, m, strand[t], \
			      (m > 1 ? st[t "_" p] : "NA"), (m > 1 ? en[t "_" p] : "NA"), \
			      st[t "_" m], en[t "_" m]
		}
	}' "${gtf}" | { read -r hdr; echo "${hdr}"; sort; } > "${ex_dir}/transcript_exons.tsv"

	echo "  transcripts: $(( $(wc -l < "${ex_dir}/transcript_exons.tsv") - 1 ))"
fi

echo "Reference data preparation completed!"
echo "Finished: $(date)"
