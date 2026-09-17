#!/bin/bash
# © 2026 GeneMolX AI LLC. All rights reserved.
# 1500 N GRANT ST STE N DENVER, CO, 80203, USA

set -euo pipefail

# GMX_clinical_variants.sh
# Clinical interpretation of germline variants from any human whole-genome sample
# Input:  curated_snps_clinical.tsv + curated_indels_clinical.tsv (48 cols)
# Output: clinical_tiers.tsv, clinical_report.txt, clinical_report.html, and
#         clinvar_plp_snps.tsv + clinvar_plp_indels.tsv (kept for compatibility)
# Depends on: GMX_annotation.sh (Step 7 must be complete)
#
# --- v2 (tiered clinical output) -------------------------------------------
# The ClinVar exact-match extraction above is retained unchanged for backward
# compatibility. v2 adds an evidence-tiered table and a human-readable report
# built from the same curated TSVs (48 cols after the Step 7 quality-field
# update: QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum).
#
# Tiers (one tier per variant, assigned by priority):
#   1 KNOWN_PATHOGENIC   ClinVar P/LP, not conflicting
#   2 RARE_LOF           rare loss-of-function (nonsense, frameshift, splice site,
#                        start-loss, stop-loss)
#   3 PREDICTED_DAMAGING rare, above ClinGen-calibrated predictor thresholds
#                        (REVEL / AlphaMissense / SpliceAI)
#   4 SECONDARY_FINDING  tier-1 variant in an ACMG SF v3.3 gene that satisfies
#                        that gene's Variants_to_Report rule
#   5 CARRIER            single heterozygous P/LP in a recessive gene
#   6 REVIEW             ClinVar conflicting classifications, or an SF-gene P/LP
#                        that fails its reporting rule - never dropped silently
#
# Outputs:
#   clinvar_plp_snps.tsv / clinvar_plp_indels.tsv  (v1, unchanged)
#   clinical_tiers.tsv                             (tiered, evidence-carrying)
#   clinical_report.txt                            (human-readable summary)
# ---------------------------------------------------------------------------

GMX_VERSION="0.1.0"
case "${1:-}" in
	-V|--version) echo "GMX_clinical_variants.sh ${GMX_VERSION} - GeneMolX WGS pipeline"; exit 0 ;;
esac

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate genemolx-wgs

echo "Started: $(date)"

# Usage: GMX_clinical_variants.sh <results_dir> [acmg_sf_table] [aligned_reads_dir] [sample_id] [--offline]
#   --offline skips the gnomAD v4.1 remote FAF lookup (uses the cache if present)
online=1
args=()
for a in "$@"; do
    if [ "$a" = "--offline" ]; then online=0; else args+=("$a"); fi
done
set -- "${args[@]}"
if [ $# -lt 1 ]; then
	echo "Usage: GMX_clinical_variants.sh <results_dir> [acmg_sf_table] [aligned_reads_dir] [sample_id] [--offline]"
	echo "  results_dir must contain curated_snps_clinical.tsv and curated_indels_clinical.tsv"
	echo "  from GMX_annotation.sh"
	exit 1
fi
[ -d "$1" ] || { echo "ERROR: results directory not found: $1"; exit 1; }

results="$(cd "$1" && pwd)"   # absolute: the gnomAD lookup runs in a subshell that changes directory
acmg_table="${2:-$(dirname "$0")/resources/acmg/ACMG_SF_v3.3_table1_full.tsv}"
aligned="${3:-$(dirname "${results}")/aligned_reads}"
sample="${4:-SAMPLE}"
gencc_table="${5:-$(dirname "$0")/resources/gencc/gencc_gene_inheritance.tsv}"
gnomad_url="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/joint"
constraint_table="${6:-$(dirname "$0")/resources/gnomad/gnomad_v4.1_constraint.tsv}"
exon_table="${7:-$(dirname "$0")/resources/gencode/transcript_exons.tsv}"

# --- thresholds ------------------------------------------------------------
RARE_MAX=0.001          # rare: population max AF <= 0.1%
BS1_AF=0.005            # BS1 for dominant / X-linked genes
BS1_AR=0.02             # BS1 where any recessive assertion exists (carrier frequencies are
                        # legitimately high; a false benign call is the costly direction)
BA1_AF=0.05             # provisional BA1
REVEL_SUP=0.644         # ClinGen calibrated: supporting
REVEL_MOD=0.773         # moderate
REVEL_STR=0.932         # strong
SPLICEAI_SUP=0.2        # ClinGen SVI: PP3 supporting
SPLICEAI_HIGH=0.5       # high-confidence splice impact
# variant-level quality review thresholds (all variants here already PASSed the
# GATK hard filters; these are the tighter thresholds used to flag calls that
# still look artefactual)
Q_DP=15; Q_GQ=30; Q_QD=5; Q_MQ=55; Q_MQRS=-2.5; Q_RPRS=-3.0
Q_FS_SNP=30; Q_FS_INDEL=20; Q_SOR=3.0
Q_VAF_HET_LO=0.25; Q_VAF_HET_HI=0.75; Q_VAF_HOM_LO=0.85

# ===========================================================================
# PART 1 - v1 output: ClinVar P/LP extraction (unchanged)
# ===========================================================================

# Columns are selected by header name (ClinVar_Significance, ClinVar_Review_Status, GT, Gene, ...)
CRV_AWK='
function getv(name) { return (name in col) ? $(col[name]) : "." }
BEGIN {
    FS = OFS = "\t"
    # Inheritance for ACMG SF genes (Gene = column 5, Inheritance = column 6)
    while ((getline line < acmg_table) > 0) {
        if (++nl == 1) continue
        split(line, a, "\t")
        if (!(a[5] in moi)) moi[a[5]] = a[6]
        else if (index(moi[a[5]], a[6]) == 0) moi[a[5]] = moi[a[5]] "/" a[6]
    }
    close(acmg_table)
}
NR == 1 {
    for (i = 1; i <= NF; i++) col[$i] = i
    print $0, "ClinVar_Stars", "BA1_Flag", "Zygosity", "Inheritance"
    next
}
{
    # Pathogenic / Likely_pathogenic, including compound values (e.g. Pathogenic_%7C_risk_factor); not Conflicting
    sig = tolower(getv("ClinVar_Significance"))
    if (sig !~ /(^|[^a-z])(likely_)?pathogenic/ || sig ~ /conflicting/) next

    # ClinVar review status -> stars
    rs = tolower(getv("ClinVar_Review_Status"))
    if (rs ~ /practice_guideline/) stars = "4"
    else if (rs ~ /reviewed_by_expert_panel/) stars = "3"
    else if (rs ~ /no_assertion|no_classification/ || rs == "" || rs == ".") stars = "0"
    else if (rs ~ /multiple_submitters/ && rs ~ /no_conflicts/) stars = "2"
    else if (rs ~ /criteria_provided/) stars = "1"
    else stars = "0"

    # BA1: gnomAD v4.1 Grpmax FAF95 > 5%
    faf = getv("gnomAD_v4.1_Grpmax_FAF95")
    ba1 = (faf != "" && faf != "." && faf + 0 > 0.05) ? "BA1" : "."

    # Zygosity
    gt = getv("GT"); gsub(/\|/, "/", gt)
    if (gt == "0/1" || gt == "1/0") zyg = "HET"
    else if (gt == "1/1") zyg = "HOM_ALT"
    else if (gt == "./." || gt == ".") zyg = "NO_CALL"
    else zyg = gt

    # Inheritance: ACMG SF table first, then Funcotator Mode_of_Inheritance
    gene = getv("Gene")
    inh = (gene in moi) ? moi[gene] : getv("Mode_of_Inheritance")
    if (inh == "") inh = "."

    print $0, stars, ba1, zyg, inh
}'

echo "Extracting Pathogenic and Likely_pathogenic SNPs..."
mawk -v acmg_table="${acmg_table}" "${CRV_AWK}" \
    "${results}/curated_snps_clinical.tsv" \
    > "${results}/clinvar_plp_snps.tsv"
snp_count=$(( $(wc -l < "${results}/clinvar_plp_snps.tsv") - 1 ))
echo "  SNPs: ${snp_count} variants"

echo "Extracting Pathogenic and Likely_pathogenic INDELs..."
mawk -v acmg_table="${acmg_table}" "${CRV_AWK}" \
    "${results}/curated_indels_clinical.tsv" \
    > "${results}/clinvar_plp_indels.tsv"
indel_count=$(( $(wc -l < "${results}/clinvar_plp_indels.tsv") - 1 ))
echo "  INDELs: ${indel_count} variants"

# QC metrics (Picard / VerifyBamID2), if the run directory is available
mcov="NA"; p20="NA"; p10="NA"; dup="NA"; freemix="NA"; aln="NA"
if [ -f "${aligned}/wgs_metrics.txt" ]; then
    read -r mcov p10 p20 <<< "$(mawk -F'\t' '/^GENOME_TERRITORY/{getline v; split(v,V,"\t"); for(i=1;i<=NF;i++) h[$i]=i; print V[h["MEAN_COVERAGE"]], V[h["PCT_10X"]], V[h["PCT_20X"]]}' "${aligned}/wgs_metrics.txt")"
fi
if [ -f "${aligned}/dedup_metrics.txt" ]; then
    dup=$(mawk -F'\t' '/^LIBRARY/{getline v; for(i=1;i<=NF;i++) h[$i]=i; split(v,V,"\t"); print V[h["PERCENT_DUPLICATION"]]}' "${aligned}/dedup_metrics.txt")
fi
if [ -f "${aligned}/verifybamid2.selfSM" ]; then
    freemix=$(mawk 'NR==2{print $7}' "${aligned}/verifybamid2.selfSM")
fi
if [ -f "${aligned}/alignment_metrics.txt" ]; then
    aln=$(mawk -F'\t' '/^CATEGORY/{for(i=1;i<=NF;i++) h[$i]=i} $1=="PAIR"{print $h["PCT_PF_READS_ALIGNED"]}' "${aligned}/alignment_metrics.txt")
fi

# -------------------
# Sample acceptance gate
# -------------------
# Germline WGS interpretation needs enough depth to genotype. Below these
# thresholds the callset is not interpretable: at low coverage the only sites
# reaching the depth filter are collapsed repeats, which produce convincing-
# looking but false loss-of-function calls. The pipeline still reports what it
# found, clearly marked, rather than silently classifying artefacts.
MIN_MEAN_COV=20         # mean genome coverage
MIN_PCT_20X=0.80        # fraction of the genome at >= 20x
MAX_FREEMIX=0.03        # contamination

qc_pass=1
qc_reasons=""
num_ok() { case "$1" in ''|NA) return 1 ;; *) return 0 ;; esac; }
if num_ok "${mcov}" && [ "$(mawk -v a="${mcov}" -v b="${MIN_MEAN_COV}" 'BEGIN{print (a<b)}')" = "1" ]; then
    qc_pass=0; qc_reasons="${qc_reasons}mean coverage ${mcov}x is below ${MIN_MEAN_COV}x; "
fi
if num_ok "${p20}" && [ "$(mawk -v a="${p20}" -v b="${MIN_PCT_20X}" 'BEGIN{print (a<b)}')" = "1" ]; then
    qc_pass=0; qc_reasons="${qc_reasons}only $(mawk -v a="${p20}" 'BEGIN{printf "%.1f", a*100}')% of the genome reaches 20x (need $(mawk -v b="${MIN_PCT_20X}" 'BEGIN{printf "%.0f", b*100}')%); "
fi
if num_ok "${freemix}" && [ "$(mawk -v a="${freemix}" -v b="${MAX_FREEMIX}" 'BEGIN{print (a>b)}')" = "1" ]; then
    qc_pass=0; qc_reasons="${qc_reasons}contamination ${freemix} exceeds ${MAX_FREEMIX}; "
fi
if [ "${qc_pass}" -eq 0 ]; then
    echo ""
    echo "############################################################"
    echo "# SAMPLE FAILED QC - RESULTS ARE NOT INTERPRETABLE"
    echo "# ${qc_reasons}"
    echo "# Variants are still listed, but none will be classified."
    echo "############################################################"
    echo ""
fi

# ===========================================================================
# PART 2 - v2 output: evidence-tiered clinical table
# ===========================================================================

echo "Building tiered clinical table..."

TIER_AWK=$(mktemp /tmp/gmx_tiers.XXXXXX.awk)
cat > "${TIER_AWK}" << 'AWKEOF'
function dec(s) {
    gsub(/%7C/, "|", s); gsub(/%2C/, ",", s); gsub(/%3B/, ";", s)
    gsub(/%3D/, "=", s); gsub(/%26/, "&", s); gsub(/%23/, "#", s)
    return s
}
function nz(s)   { return (s == "" || s == ".") ? "" : s }
function num(s)  { return (nz(s) == "") ? -1 : s + 0 }
function getv(n) { return (n in col) ? $(col[n]) : "" }
function put(v)  { return (nz(v) == "") ? "." : dec(v) }
function par(p)  { return ((p >= 10001 && p <= 2781479) || (p >= 155701383 && p <= 156030895)) }
function tname(t) {
    if (t == 1) return "KNOWN_PATHOGENIC"
    if (t == 2) return "RARE_LOF"
    if (t == 3) return "PREDICTED_DAMAGING"
    if (t == 4) return "SECONDARY_FINDING"
    if (t == 5) return "CARRIER"
    if (t == 6) return "REVIEW"
    return "NA"
}
BEGIN {
    FS = OFS = "\t"
    # ACMG SF v3.3: Gene=5, Inheritance=6, Variants_to_Report=7, Phenotype=2, MIM=4
    while ((getline line < acmg_table) > 0) {
        if (++nl == 1) continue
        split(line, a, "\t"); g = a[5]
        if (!(g in sf_rule)) { sf_rule[g] = a[7]; sf_inh[g] = a[6]; sf_dis[g] = a[2] }
        else {
            if (index(sf_inh[g], a[6]) == 0) sf_inh[g] = sf_inh[g] "/" a[6]
            if (index(sf_dis[g], a[2]) == 0) sf_dis[g] = sf_dis[g] "; " a[2]
        }
    }
    close(acmg_table)
    # GenCC (CC0): Gene=1, Inheritance=2, Validity=3, Disease=4
    if (gencc_table != "") {
        nl = 0
        while ((getline line < gencc_table) > 0) {
            if (++nl == 1) continue
            split(line, a, "\t")
            gc_inh[a[1]] = a[2]; gc_val[a[1]] = a[3]; gc_dis[a[1]] = a[4]
        }
        close(gencc_table)
    }
    n = split("NONSENSE FRAME_SHIFT_INS FRAME_SHIFT_DEL SPLICE_SITE START_CODON_SNP START_CODON_DEL START_CODON_INS NONSTOP", L, " ")
    for (i = 1; i <= n; i++) lof[L[i]] = 1
    nr = 0; nx = 0; hx = 0; ny = 0
}
FNR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
{
    # --- genotype (needed for the sex check on every record) --------------
    gt = getv("GT"); zgt = gt; gsub(/\|/, "/", zgt)
    if (zgt == "0/1" || zgt == "1/0") zyg = "HET"
    else if (zgt == "1/1") zyg = "HOM_ALT"
    else if (zgt == "./." || zgt == ".") zyg = "NO_CALL"
    else zyg = zgt
    if ($1 == "chrX" && !par($2)) { nx++; if (zyg == "HET") hx++ }
    else if ($1 == "chrY" && !par($2)) ny++

    gene = nz(getv("Gene")); if (gene == "") gene = "."
    cls  = getv("Variant_Classification")
    vtyp = getv("Variant_Type")

    # --- ClinVar ----------------------------------------------------------
    sig  = dec(getv("ClinVar_Significance")); lsig = tolower(sig)
    isplp  = (lsig ~ /(^|[^a-z])(likely_)?pathogenic/ && lsig !~ /conflicting/)
    isconf = (lsig ~ /conflicting/)
    isben  = (!isplp && !isconf && lsig ~ /benign/)
    rs = tolower(getv("ClinVar_Review_Status"))
    if (rs ~ /practice_guideline/) st = 4
    else if (rs ~ /reviewed_by_expert_panel/) st = 3
    else if (rs ~ /no_assertion|no_classification/ || rs == "" || rs == ".") st = 0
    else if (rs ~ /multiple_submitters/ && rs ~ /no_conflicts/) st = 2
    else if (rs ~ /criteria_provided/) st = 1
    else st = 0

    # --- population frequency (gnomAD v4.1 if present, else v2.1) ----------
    f41 = num(getv("gnomAD_v4.1_Grpmax_FAF95"))
    pe  = num(getv("gnomAD_Exome_AF_PopMax")); pg = num(getv("gnomAD_Genome_AF_PopMax"))
    ae  = num(getv("gnomAD_Exome_AF"));        ag = num(getv("gnomAD_Genome_AF"))
    if (f41 >= 0)                { freq = f41; fsrc = "gnomADv4.1_FAF95" }
    else if (pe >= 0 || pg >= 0) { freq = (pe > pg ? pe : pg); fsrc = "gnomADv2.1_popmax" }
    else if (ae >= 0 || ag >= 0) { freq = (ae > ag ? ae : ag); fsrc = "gnomADv2.1_AF" }
    else                         { freq = -1;  fsrc = "absent" }
    rare = (freq < 0 || freq <= rare_max)
    if (freq > ba1_af)      fflag = "BA1_provisional"
    else if (freq > bs1_af) fflag = "BS1_provisional"
    else if (freq < 0)      fflag = "PM2_absent"
    else                    fflag = "."

    # --- in silico --------------------------------------------------------
    rev = num(getv("REVEL_Score")); am = tolower(dec(getv("AlphaMissense_Class")))
    spl = num(getv("SpliceAI_DS_Max"))
    ev = ""
    if (rev >= revel_str)      ev = ev ";REVEL=" rev "(strong)"
    else if (rev >= revel_mod) ev = ev ";REVEL=" rev "(moderate)"
    else if (rev >= revel_sup) ev = ev ";REVEL=" rev "(supporting)"
    if (am ~ /likely_pathogenic/) ev = ev ";AlphaMissense=likely_pathogenic"
    if (spl >= splice_high)     ev = ev ";SpliceAI=" spl "(high)"
    else if (spl >= splice_sup) ev = ev ";SpliceAI=" spl "(supporting)"
    pred = (rev >= revel_sup || am ~ /likely_pathogenic/ || spl >= splice_sup)

    vaf = num(getv("VAF")); dp = num(getv("DP")); gq = num(getv("GQ"))

    # --- variant-level quality flags --------------------------------------
    qd = num(getv("QD")); fs = num(getv("FS")); sor = num(getv("SOR"))
    mq = num(getv("MQ")); mqrs = getv("MQRankSum"); rprs = getv("ReadPosRankSum")
    q = ""; nq = 0
    if (dp >= 0 && dp < q_dp)        { q = q ";DP=" dp;   nq++ }
    if (gq >= 0 && gq < q_gq)        { q = q ";GQ=" gq;   nq++ }
    if (qd >= 0 && qd < q_qd)        { q = q ";QD=" qd;   nq++ }
    if (mq >= 0 && mq < q_mq)        { q = q ";MQ=" mq;   nq++ }
    if (nz(mqrs) != "" && mqrs + 0 < q_mqrs) { q = q ";MQRankSum=" mqrs; nq++ }
    if (nz(rprs) != "" && rprs + 0 < q_rprs) { q = q ";ReadPosRankSum=" rprs; nq++ }
    fs_lim = (vtyp == "SNP") ? q_fs_snp : q_fs_indel
    if (fs >= 0 && fs > fs_lim)      { q = q ";FS=" fs;   nq++ }
    if (sor >= 0 && sor > q_sor)     { q = q ";SOR=" sor; nq++ }
    if (zyg == "HET" && vaf >= 0 && (vaf < q_vaf_het_lo || vaf > q_vaf_het_hi)) { q = q ";VAF=" vaf; nq++ }
    if (zyg == "HOM_ALT" && vaf >= 0 && vaf < q_vaf_hom_lo)                     { q = q ";VAF=" vaf; nq++ }
    if ($1 == "chrM" || $1 == "chrMT") { q = q ";chrM_diploid_caller"; nq++ }
    # depth far above the genome mean indicates a collapsed repeat, where reads
    # from several genomic copies stack at one locus
    if (mean_cov > 0 && dp > 0 && dp > 5 * mean_cov) { q = q ";DP_outlier=" dp "x_vs_" mean_cov "x_mean"; nq++ }

    # --- tier -------------------------------------------------------------
    if (isplp)                       { tier = 1; why = "ClinVar " sig " (" st " star)" }
    else if (isconf)                 { tier = 6; why = "ClinVar conflicting classifications" }
    else if (rare && !isben && (cls in lof)) { tier = 2; why = "rare " cls }
    else if (rare && !isben && pred) { tier = 3; why = "rare " cls ", predictor support" }
    else next

    sub(/^;/, "", ev); if (ev == "") ev = "."
    # inheritance: ACMG SF table -> GenCC -> Funcotator Mode_of_Inheritance
    if (gene in sf_inh)        inh = sf_inh[gene]
    else if (gene in gc_inh && gc_inh[gene] != ".") inh = gc_inh[gene]
    else                       inh = nz(getv("Mode_of_Inheritance"))
    if (inh == "") inh = "."
    validity = (gene in gc_val) ? gc_val[gene] : "."
    acmg_sf = (gene in sf_rule) ? "True" : "False"
    if (gene in sf_dis) acmg_dis = sf_dis[gene]
    else if (put(getv("ACMG_Disease")) != ".") acmg_dis = put(getv("ACMG_Disease"))
    else if (gene in gc_dis) acmg_dis = gc_dis[gene]
    else acmg_dis = "."

    nr++
    if (tier == 2) lofgene[gene]++
    t_tier[nr] = tier; t_why[nr] = why; t_gene[nr] = gene; t_zyg[nr] = zyg
    t_plp[nr] = isplp; t_lof[nr] = (cls in lof) ? 1 : 0
    t_prot[nr] = put(getv("Protein_Change")); t_inh[nr] = inh
    t_chr[nr] = $1; t_q[nr] = q; t_nq[nr] = nq; t_val[nr] = validity
    t_row[nr] = gene OFS $1 OFS $2 OFS $4 OFS $5 OFS "@ZYG@" OFS cls OFS t_prot[nr] \
        OFS put(getv("cDNA_Change")) OFS put(getv("Transcript")) OFS put(getv("Exon")) \
        OFS "@WHY@" OFS ev OFS (nz(sig) == "" ? "." : sig) OFS st OFS put(getv("ClinVar_Disease")) \
        OFS acmg_sf OFS "@RULE@" OFS acmg_dis OFS inh OFS validity \
        OFS (freq < 0 ? "NA" : freq) OFS fsrc OFS fflag \
        OFS (rev < 0 ? "." : rev) OFS put(getv("AlphaMissense_Class")) OFS (spl < 0 ? "." : spl) \
        OFS "@CONF@" OFS "@QC@" \
        OFS gt OFS put(getv("DP")) OFS put(getv("AD")) OFS put(getv("VAF")) OFS put(getv("GQ")) \
        OFS put(getv("QUAL")) OFS put(getv("QD")) OFS put(getv("FS")) OFS put(getv("SOR")) \
        OFS put(getv("MQ")) OFS put(getv("MQRankSum")) OFS put(getv("ReadPosRankSum")) \
        OFS put(getv("dbSNP_RS")) OFS put(getv("Genome_Change")) OFS put(getv("ClinVar_HGVS"))
}
END {
    # --- genetic sex from chrX heterozygosity (non-PAR) --------------------
    xhet = (nx > 0) ? hx / nx : -1
    if (xhet < 0)         sex = "UNKNOWN"
    else if (xhet < 0.15) sex = "XY"
    else if (xhet > 0.30) sex = "XX"
    else                  sex = "UNCERTAIN"
    printf "#SEX\t%s\t%.4f\t%d\t%d\n", sex, (xhet < 0 ? 0 : xhet), nx, ny > "/dev/stderr"

    # allele count per gene, to apply recessive / two-variant reporting rules
    for (i = 1; i <= nr; i++)
        if (t_plp[i]) plp[t_gene[i]] += (t_zyg[i] == "HOM_ALT" ? 2 : 1)

    for (i = 1; i <= nr; i++) {
        t = t_tier[i]; g = t_gene[i]; why = t_why[i]; rule = "."
        zyg = t_zyg[i]; q = t_q[i]; nq = t_nq[i]

        # sex-aware zygosity and sex-inconsistent call flags
        if (sex == "XY" && (t_chr[i] == "chrX" || t_chr[i] == "chrY")) {
            if (zyg == "HOM_ALT") zyg = "HEMI"
            else if (zyg == "HET") { q = q ";het_on_" t_chr[i] "_in_XY_sample"; nq++ }
        }
        if (sex == "XX" && t_chr[i] == "chrY") { q = q ";chrY_call_in_XX_sample"; nq++ }

        # loss-of-function pile-up in one gene: assembly / mapping artefact pattern
        if (t_tier[i] == 2 && lofgene[g] >= 3) { q = q ";LoF_cluster=" lofgene[g]; nq++ }

        sub(/^;/, "", q); if (q == "") q = "."
        conf = (nq == 0) ? "HIGH" : ((nq == 1) ? "MEDIUM" : "LOW")

        if (t == 1) {
            if (g in sf_rule) {
                rule = sf_rule[g]; ok = 1; note = ""
                if (rule ~ /2 variants|two het/ && plp[g] < 2) {
                    ok = 0; note = "single P/LP allele in recessive SF gene - not a reportable secondary finding"
                }
                if (rule ~ /C282Y/ && !(zyg == "HOM_ALT" && t_prot[i] ~ /C282Y/)) {
                    ok = 0; note = "HFE reportable only as p.C282Y homozygote"
                }
                if (rule ~ /truncating only/ && !t_lof[i]) {
                    ok = 0; note = "SF gene reportable for truncating variants only"
                }
                if (ok) { t = 4; why = why "; ACMG SF v3.3 reportable" }
                else if ((t_inh[i] ~ /AR/ && zyg == "HET") || (t_inh[i] ~ /XL/ && zyg == "HET" && sex == "XX")) { t = 5; why = why "; " note }
                else { t = 6; why = why "; " note }
            } else if (t_inh[i] ~ /AR/ && zyg == "HET" && plp[g] < 2) {
                t = 5; why = why "; single heterozygote in recessive gene"
            } else if (t_inh[i] ~ /XL/ && zyg == "HET" && sex == "XX") {
                t = 5; why = why "; heterozygous in X-linked gene (XX sample)"
            } else if (t_inh[i] == "." && zyg == "HET") {
                why = why "; inheritance unknown - carrier status not assessable"
            }
        }
        line = t_row[i]
        sub(/@WHY@/, why, line); sub(/@RULE@/, rule, line)
        sub(/@ZYG@/, zyg, line); sub(/@CONF@/, conf, line); sub(/@QC@/, q, line)
        crank = (conf == "HIGH") ? 1 : ((conf == "MEDIUM") ? 2 : 3)
        # unnamed Ensembl gene models rank below named genes - no disease association
        if (g ~ /^ENSG[0-9]/) crank += 5
        # genes without established gene-disease validity rank below those with it
        if (t_val[i] !~ /Definitive|Strong|Moderate/) crank += 2
        if (sex == "XX" && t_chr[i] == "chrY") crank += 3
        print crank, t, tname(t), line
    }
}
AWKEOF

TIER_HEADER="Tier\tTier_Name\tGene\tCHROM\tPOS\tREF\tALT\tZygosity\tVariant_Classification\tProtein_Change\tcDNA_Change\tTranscript\tExon\tReason\tEvidence\tClinVar_Significance\tClinVar_Stars\tClinVar_Disease\tACMG_SF\tACMG_Report_Rule\tACMG_Disease\tInheritance\tGene_Validity\tFreq_Max\tFreq_Source\tFreq_Flag\tREVEL_Score\tAlphaMissense_Class\tSpliceAI_DS_Max\tConfidence\tQC_Flags\tGT\tDP\tAD\tVAF\tGQ\tQUAL\tQD\tFS\tSOR\tMQ\tMQRankSum\tReadPosRankSum\tdbSNP_RS\tGenome_Change\tClinVar_HGVS"

printf "${TIER_HEADER}\n" > "${results}/clinical_tiers.tsv"
mawk -v acmg_table="${acmg_table}" -v gencc_table="${gencc_table}" \
     -v rare_max="${RARE_MAX}" -v bs1_af="${BS1_AF}" -v ba1_af="${BA1_AF}" \
     -v revel_sup="${REVEL_SUP}" -v revel_mod="${REVEL_MOD}" -v revel_str="${REVEL_STR}" \
     -v splice_sup="${SPLICEAI_SUP}" -v splice_high="${SPLICEAI_HIGH}" \
     -v q_dp="${Q_DP}" -v q_gq="${Q_GQ}" -v q_qd="${Q_QD}" -v q_mq="${Q_MQ}" \
     -v q_mqrs="${Q_MQRS}" -v q_rprs="${Q_RPRS}" -v q_fs_snp="${Q_FS_SNP}" \
     -v q_fs_indel="${Q_FS_INDEL}" -v q_sor="${Q_SOR}" \
     -v q_vaf_het_lo="${Q_VAF_HET_LO}" -v q_vaf_het_hi="${Q_VAF_HET_HI}" -v q_vaf_hom_lo="${Q_VAF_HOM_LO}" \
     -v mean_cov="$(num_ok "${mcov}" && echo "${mcov}" || echo 0)" \
     -f "${TIER_AWK}" \
     "${results}/curated_snps_clinical.tsv" "${results}/curated_indels_clinical.tsv" \
     2> "${results}/.sexcheck" \
  | sort -t$'\t' -k2,2n -k1,1n -k4,4 -k5,5 -k6,6n \
  | cut -f2- >> "${results}/clinical_tiers.tsv"

rm -f "${TIER_AWK}"

tier_total=$(( $(wc -l < "${results}/clinical_tiers.tsv") - 1 ))
echo "  Tiered rows: ${tier_total}"

# ===========================================================================
# PART 2.5 - gnomAD v4.1 filtering allele frequency (BA1 / BS1)
# ===========================================================================
# The gnomAD v4.1 joint release is ~1 TB, so it is not mirrored locally. Only
# the tiered variants (a few hundred) are queried, by remote tabix, and cached.
# grpmax FAF95 (fafmax_faf95_max_joint) is the frequency statistic ACMG/AMP
# BA1 and BS1 are defined on.

faf_cache="${results}/gnomad_v41_faf.tsv"
[ -f "${faf_cache}" ] || printf "KEY\tAF_joint\tFAF95\n" > "${faf_cache}"

# remote index files are cached here, not in the caller's working directory
tbi_cache="${results}/.tabix_cache"
mkdir -p "${tbi_cache}"

if [ "${online}" -eq 1 ]; then
    echo "Querying gnomAD v4.1 (joint) for tiered variants..."
    reg="${results}/.faf_regions.bed"
    mawk -F'\t' -v cache="${faf_cache}" '
        BEGIN { while ((getline l < cache) > 0) { split(l, a, "\t"); have[a[1]] = 1 } }
        NR > 1 { k = $4 ":" $5 ":" $6 ":" $7; qk = $4 ":" $5 ":QUERIED:."
                 if (!(k in have) && !(qk in have) && !(k in seen)) { seen[k] = 1; print $4, $5 - 1, $5 } }
    ' OFS='\t' "${results}/clinical_tiers.tsv" | sort -u -k1,1 -k2,2n > "${reg}"
    need=$(wc -l < "${reg}")
    echo "  variants needing lookup: ${need}"
    if [ "${need}" -gt 0 ]; then
        # remote reads drop on long streams: query in small chunks and retry
        parse_faf='{
            af = "."; faf = "."
            n = split($8, I, ";")
            for (i = 1; i <= n; i++) {
                if (I[i] ~ /^AF_joint=/)                    { sub(/^AF_joint=/, "", I[i]);                 af = I[i] }
                else if (I[i] ~ /^fafmax_faf95_max_joint=/) { sub(/^fafmax_faf95_max_joint=/, "", I[i]);   faf = I[i] }
            }
            print $1 ":" $2 ":" $4 ":" $5, af, faf
        }'
        for c in $(cut -f1 "${reg}" | sort -u); do
            mawk -v c="${c}" '$1 == c' "${reg}" > "${results}/.faf_chrom.bed"
            split -l 10 "${results}/.faf_chrom.bed" "${results}/.faf_chunk_"
            failed=0
            for chunk in "${results}"/.faf_chunk_*; do
                ok=0
                for try in 1 2 3 4 5; do
                    if ( cd "${tbi_cache}" && timeout 600 tabix "${gnomad_url}/gnomad.joint.v4.1.sites.${c}.vcf.bgz" \
                           -R "${chunk}" 2> /dev/null ) > "${results}/.faf_raw"; then
                        mawk -F'\t' "${parse_faf}" OFS='\t' "${results}/.faf_raw" >> "${faf_cache}"
                        mawk '{ print $1 ":" $3 ":QUERIED:.", ".", "." }' OFS='\t' "${chunk}" >> "${faf_cache}"
                        ok=1; break
                    fi
                    sleep 3
                done
                [ "${ok}" -eq 1 ] || failed=$(( failed + 1 ))
            done
            rm -f "${results}"/.faf_chunk_*
            [ "${failed}" -eq 0 ] || echo "  warning: ${failed} chunk(s) failed for ${c} (local frequencies kept there)"
        done
        rm -f "${results}/.faf_raw"
        # de-duplicate the cache (retries and partial streams can repeat keys)
        { head -1 "${faf_cache}"; tail -n +2 "${faf_cache}" | sort -u; } > "${faf_cache}.tmp"
        mv -f "${faf_cache}.tmp" "${faf_cache}"
    fi
    rm -f "${reg}" "${results}/.faf_chrom.bed"
fi

# apply the cached FAF95 values: update frequency columns, and move variants that
# turn out to be common in gnomAD v4.1 out of the candidate tiers into review
mawk -F'\t' -v cache="${faf_cache}" -v rare_max="${RARE_MAX}" -v bs1_af="${BS1_AF}" -v ba1_af="${BA1_AF}" -v ranked="${results}/.tiers_ranked" '
BEGIN {
    OFS = "\t"
    while ((getline l < cache) > 0) { split(l, a, "\t"); if (a[1] != "KEY") { gaf[a[1]] = a[2]; gfaf[a[1]] = a[3] } }
    close(cache)
}
NR == 1 { print; next }
{
    k = $4 ":" $5 ":" $6 ":" $7
    if (k in gfaf) {
        f = (gfaf[k] != "." && gfaf[k] != "") ? gfaf[k] + 0 : ((gaf[k] != "." && gaf[k] != "") ? gaf[k] + 0 : -1)
        if (f >= 0) {
            $24 = f
            $25 = (gfaf[k] != "." && gfaf[k] != "") ? "gnomADv4.1_FAF95" : "gnomADv4.1_AF_joint"
            if (f > ba1_af)      $26 = "BA1"
            else if (f > bs1_af) $26 = "BS1"
            else                 $26 = "."
            if (($1 == 2 || $1 == 3) && f > rare_max) {
                $14 = $14 "; not rare in gnomAD v4.1 (FAF95=" f ")"
                $1 = 6; $2 = "REVIEW"
            }
        }
    } else if ($25 ~ /^gnomADv2/) {
        $26 = ($26 == ".") ? "." : $26   # unchanged; v2.1 flags stay provisional
    }
    # re-derive the sort rank after any tier change
    crank = ($30 == "HIGH") ? 1 : (($30 == "MEDIUM") ? 2 : 3)
    if ($3 ~ /^ENSG[0-9]/) crank += 5
    if ($23 !~ /Definitive|Strong|Moderate/) crank += 2
    print crank, $0 > ranked
}' "${results}/clinical_tiers.tsv" > "${results}/.tiers_head"

{ cat "${results}/.tiers_head"
  sort -t$'\t' -k2,2n -k1,1n -k4,4 -k5,5 -k6,6n "${results}/.tiers_ranked" | cut -f2- ; } > "${results}/clinical_tiers.tsv"
rm -f "${results}/.tiers_ranked" "${results}/.tiers_head"

faf_n=$(( $(wc -l < "${faf_cache}") - 1 ))
echo "  gnomAD v4.1 records cached: ${faf_n}"

# ===========================================================================
# PART 2.6 - ACMG/AMP classification (point system)
# ===========================================================================
# Evidence codes are assigned from the data actually available in this pipeline
# and combined with the Bayesian point scale used by the ClinGen Sequence
# Variant Interpretation working group (Tavtigian 2020):
#   very strong 8, strong 4, moderate 2, supporting 1; benign codes negative.
#   Pathogenic >= 10, Likely pathogenic 6-9, VUS 0-5, Likely benign -1..-6,
#   Benign <= -7, and BA1 is stand-alone benign.
#
# Implemented: PVS1 (with NMD-escape downgrade), PM2_Supporting, PM4, PP2,
#              PP3 / BP4 (ClinGen-calibrated REVEL and SpliceAI), BA1, BS1, BP7.
# Not implemented (no data): PS1, PS2/PM6, PS3/BS3, PS4, PM1, PM3, PM5, PP1/BS4,
#              BS2, BP1, BP2, BP3, BP5. PP5 and BP6 are retired by ClinGen.

echo "Classifying variants (ACMG/AMP point system)..."

CLASS_AWK=$(mktemp /tmp/gmx_acmg.XXXXXX.awk)
cat > "${CLASS_AWK}" << 'AWKEOF'
function add(code, points) { codes = (codes == "" ? code : codes ";" code); pts += points }
function n(x) { return (x == "" || x == "." || x == "NA") ? -1 : x + 0 }
BEGIN {
    FS = OFS = "\t"
    while ((getline line < constraint_table) > 0) {
        if (++nl == 1) continue
        split(line, a, "\t"); loeuf[a[1]] = a[2]; misz[a[1]] = a[3]; pli[a[1]] = a[4]
    }
    close(constraint_table)
    nl = 0
    while ((getline line < exon_table) > 0) {
        if (++nl == 1) continue
        split(line, a, "\t")
        ex_n[a[1]] = a[2]; ex_strand[a[1]] = a[3]
        ex_ps[a[1]] = a[4]; ex_pe[a[1]] = a[5]; ex_ls[a[1]] = a[6]; ex_le[a[1]] = a[7]
    }
    close(exon_table)
    split("NONSENSE FRAME_SHIFT_INS FRAME_SHIFT_DEL SPLICE_SITE START_CODON_SNP START_CODON_DEL START_CODON_INS", L, " ")
    for (i in L) null_var[L[i]] = 1
}
NR == 1 { print $0, "ACMG_Class", "ACMG_Points", "ACMG_Codes", "Gene_LOEUF", "Gene_MisZ"; next }
{
    gene = $3; cls = $9; tx = $12; sub(/\..*$/, "", tx)
    pos = $5 + 0; exon = n($13)
    freq = ($24 == "NA") ? -1 : $24 + 0; fsrc = $25
    revel = n($27); spl = n($29); validity = $23; inh = $22
    lo = (gene in loeuf) ? n(loeuf[gene]) : -1
    mz = (gene in misz)  ? n(misz[gene])  : -1
    pl = (gene in pli)   ? n(pli[gene])   : -1
    codes = ""; pts = 0; standalone = ""

    # --- population frequency ---------------------------------------------
    # only gnomAD v4.1 grpmax FAF95 is used for BA1/BS1; a v2.1 fallback is not
    # strong enough to carry benign evidence
    v41 = (fsrc ~ /v4\.1/)
    if (v41 && freq > ba1_af)       { standalone = "BA1"; add("BA1", 0) }
    else if (v41 && freq > ((inh ~ /AR/) ? bs1_ar : bs1_af)) add("BS1", -4)
    else if ((freq < 0 && fsrc == "absent") || (v41 && freq >= 0 && freq < 1e-5)) add("PM2_Supporting", 1)

    # --- PVS1: null variant in a gene where loss of function causes disease -
    if (cls in null_var) {
        # constraint reflects selection against heterozygous LoF, so it cannot
        # gate PVS1 in recessive genes
        lof_gene = (validity ~ /Definitive|Strong/) && (inh ~ /AR/ || lo >= 0 && lo < 0.6 || pl >= 0.9)
        if (lof_gene) {
            nmd_escape = 0
            if (tx in ex_n && exon > 0) {
                if (exon >= ex_n[tx] + 0) nmd_escape = 1
                else if (exon == ex_n[tx] - 1 && ex_ps[tx] != "NA") {
                    d = (ex_strand[tx] == "+") ? (ex_pe[tx] + 0 - pos) : (pos - (ex_ps[tx] + 0))
                    if (d >= 0 && d <= 50) nmd_escape = 1
                }
            }
            if (cls ~ /START_CODON/)  add("PVS1_Moderate", 2)
            else if (nmd_escape)      add("PVS1_Moderate", 2)
            else if (tx in ex_n)      add("PVS1", 8)
            else                      add("PVS1_Strong", 4)      # NMD status unknown
        } else {
            codes = (codes == "" ? "" : codes ";") "PVS1_not_applied(LoF_mechanism_not_established)"
        }
    }

    # --- PM4: protein length change ---------------------------------------
    if (cls ~ /IN_FRAME/ || cls == "NONSTOP") add("PM4", 2)

    # --- PP2: missense in a missense-constrained disease gene --------------
    if (cls == "MISSENSE" && mz >= 3.09 && validity ~ /Definitive|Strong/) add("PP2", 1)

    # --- PP3 / BP4: one computational code only, strongest available -------
    p3 = 0; p3c = ""
    if (cls == "MISSENSE" && revel >= 0) {
        if (revel >= 0.932)      { p3 = 4; p3c = "PP3_Strong" }
        else if (revel >= 0.773) { p3 = 2; p3c = "PP3_Moderate" }
        else if (revel >= 0.644) { p3 = 1; p3c = "PP3_Supporting" }
        else if (revel <= 0.183) { p3 = -2; p3c = "BP4_Moderate" }
        else if (revel <= 0.290) { p3 = -1; p3c = "BP4_Supporting" }
    }
    if (spl >= 0.5 && p3 < 2)       { p3 = 2; p3c = "PP3_Moderate(splice)" }
    else if (spl >= 0.2 && p3 < 1)  { p3 = 1; p3c = "PP3_Supporting(splice)" }
    if (p3c != "") add(p3c, p3)

    # --- BP7: silent or intronic with no predicted splice effect ----------
    # BP7 also requires the nucleotide to be non-conserved, which is not
    # annotated here, so it is applied to synonymous variants only
    if (cls == "SILENT" && spl >= 0 && spl <= 0.1) add("BP7", -1)

    # --- combine -----------------------------------------------------------
    if (standalone == "BA1")   klass = "Benign"
    else if (pts >= 10)        klass = "Pathogenic"
    else if (pts >= 6)         klass = "Likely_pathogenic"
    else if (pts >= 0)         klass = "VUS"
    else if (pts >= -6)        klass = "Likely_benign"
    else                       klass = "Benign"
    # a variant that fails quality review is not classified; laboratories
    # confirm or exclude such a call before interpreting it
    if (qc_pass == 0)      klass = "Not_classified(sample_failed_QC)"
    else if ($30 == "LOW") klass = "Not_classified(low_quality_call)"
    if (codes == "") codes = "."
    print $0, klass, pts, codes, (lo < 0 ? "." : lo), (mz < 0 ? "." : mz)
}
AWKEOF

mawk -v constraint_table="${constraint_table}" -v exon_table="${exon_table}" \
     -v ba1_af="${BA1_AF}" -v bs1_af="${BS1_AF}" -v bs1_ar="${BS1_AR}" -v qc_pass="${qc_pass}" \
     -f "${CLASS_AWK}" "${results}/clinical_tiers.tsv" > "${results}/.tiers_classified"
mv -f "${results}/.tiers_classified" "${results}/clinical_tiers.tsv"
rm -f "${CLASS_AWK}"

echo "  $(mawk -F'\t' 'NR>1{c[$47]++} END{for(k in c) printf "%s=%d ", k, c[k]}' "${results}/clinical_tiers.tsv")"

# ===========================================================================
# PART 3 - human-readable report
# ===========================================================================

echo "Writing clinical report..."
REPORT="${results}/clinical_report.txt"

cnt() { mawk -F'\t' -v t="$1" 'NR>1 && $1==t' "${results}/clinical_tiers.tsv" | wc -l; }
# gene | zygosity | HGVS | ClinVar (stars) | frequency + flag | confidence | reason
show() {
    mawk -F'\t' -v t="$1" -v lim="${2:-0}" -v sf_first="${3:-0}" '
    NR > 1 && $1 == t {
        n++
        if (sf_first && $19 != "True" && pass == 1) { hold[++nh] = $0; next }
        if (lim > 0 && shown >= lim) { more++; next }
        emit($0); shown++
    }
    function reason(r) { sub(/^ClinVar [^;]*; */, "", r); return substr(r, 1, 58) }
    function emit(L,   a, h, f) {
        split(L, a, "\t")
        h = (a[10] != "." ? a[10] : a[11])
        f = (a[24] == "NA" ? "NA" : substr(a[24], 1, 9))
        if (a[26] != "." && a[26] != "PM2_absent") f = f " " a[26]
        printf "  %-10s %-7s %-20s %-26s %-24s %-18s %-6s %s\n", \
            a[3], a[8], substr(h, 1, 20), substr(a[47], 1, 24) " (" a[48] ")", \
            substr(a[16], 1, 22) " (" a[17] "*)", f, a[30], substr(reason(a[14]), 1, 34)
    }
    END {
        for (i = 1; i <= nh; i++) {
            if (lim > 0 && shown >= lim) { more++; continue }
            emit(hold[i]); shown++
        }
        if (more) printf "  ... and %d more (see clinical_tiers.tsv)\n", more
        if (!n) print "  none"
    }' pass=1 "${results}/clinical_tiers.tsv"
}

sexline="not determined"
if [ -f "${results}/.sexcheck" ]; then
    sexline=$(mawk -F'\t' '$1=="#SEX"{printf "%s (chrX non-PAR het %.1f%% of %d calls; chrY calls %d)", $2, $3*100, $4, $5}' "${results}/.sexcheck")
    rm -f "${results}/.sexcheck"
fi
plp_no_inh=$(mawk -F'\t' 'NR>1 && ($1==1) && $22=="." && $8=="HET"' "${results}/clinical_tiers.tsv" | wc -l)

{
echo "==========================================================================="
echo " GeneMolX - Clinical Variant Report"
echo " Sample: ${sample}"
echo " Generated: $(date)"
echo " Pipeline: GMX_variant_calling.sh -> GMX_annotation.sh -> GMX_clinical_variants.sh"
echo "==========================================================================="
echo ""
if [ "${qc_pass}" -eq 0 ]; then
echo "***************************************************************************"
echo " SAMPLE FAILED QUALITY CONTROL - THESE RESULTS ARE NOT INTERPRETABLE"
echo " ${qc_reasons}"
echo " No variant has been classified. At this depth the only variants passing"
echo " the depth filter are pile-ups in collapsed repeats, which look like"
echo " convincing loss-of-function calls but are artefacts. Re-sequence to the"
echo " required depth before drawing any conclusion from this report."
echo "***************************************************************************"
echo ""
fi
echo "1. RESULT SUMMARY"
echo "-----------------"
printf "   Tier 4  Secondary findings (ACMG SF v3.3 reportable) : %s\n" "$(cnt 4)"
printf "   Tier 1  Known pathogenic (ClinVar P/LP)              : %s\n" "$(cnt 1)"
printf "   Tier 5  Carrier (single het in recessive gene)       : %s\n" "$(cnt 5)"
printf "   Tier 2  Rare loss-of-function candidates             : %s\n" "$(cnt 2)"
printf "   Tier 3  Rare predicted-damaging candidates           : %s\n" "$(cnt 3)"
printf "   Tier 6  Flagged for review (conflicting / rule-fail) : %s\n" "$(cnt 6)"
echo ""
printf "   ACMG/AMP classification of the %s tiered variants:\n" "${tier_total}"
mawk -F'\t' 'NR>1 { k = $47; sub(/\(low_quality_call\)/, "", k); c[k]++; if ($47 ~ /low_quality/) lq++ }
     END { split("Pathogenic Likely_pathogenic VUS Likely_benign Benign", o, " ")
           for (i = 1; i <= 5; i++) if (o[i] in c) printf "     %-20s %d\n", o[i], c[o[i]]
           if (lq) printf "     (of which %d are low-quality calls)\n", lq }' "${results}/clinical_tiers.tsv"
echo ""
echo "2. QUALITY CONTROL"
echo "------------------"
printf "   Mean coverage        : %s\n" "${mcov}"
printf "   Genome >=10x / >=20x : %s / %s\n" "${p10}" "${p20}"
printf "   Duplication rate     : %s\n" "${dup}"
printf "   Contamination (FREEMIX): %s\n" "${freemix}"
printf "   Reads aligned        : %s\n" "${aln}"
echo ""
printf "   Genetic sex check    : %s\n" "${sexline}"
echo ""
echo "3. TIER 4 - SECONDARY FINDINGS (ACMG SF v3.3)"
echo "---------------------------------------------"
show 4
echo ""
echo "4. TIER 1 - KNOWN PATHOGENIC"
echo "----------------------------"
show 1
echo ""
echo "5. TIER 5 - CARRIER STATUS"
echo "--------------------------"
show 5
if [ "${plp_no_inh}" -gt 0 ]; then
    echo "   NOTE: ${plp_no_inh} heterozygous P/LP variant(s) have no mode of inheritance in"
    echo "   ACMG SF v3.3 or GenCC, so carrier status could not be assigned. It is NOT"
    echo "   excluded for them - the gene has no curated gene-disease assertion."
fi
echo ""
echo "6. TIER 2 - RARE LOSS-OF-FUNCTION CANDIDATES (top 25 shown)"
echo "-----------------------------------------------------------"
show 2 25
echo ""
echo "7. TIER 3 - RARE PREDICTED-DAMAGING CANDIDATES (top 25 shown)"
echo "-------------------------------------------------------------"
show 3 25
echo ""
echo "8. TIER 6 - FLAGGED FOR REVIEW (ACMG SF genes listed first, top 25)"
echo "-------------------------------------------------------------------"
show 6 25 1
echo ""
echo "9. ACMG/AMP vs ClinVar - DISCORDANT VARIANTS"
echo "--------------------------------------------"
echo "   Rows where the automated classification and the ClinVar submission point"
echo "   in different directions. Both are shown; neither overrides the other."
mawk -F'\t' '
NR > 1 {
    cv = tolower($16); k = $47
    plp = (cv ~ /(^|[^a-z])(likely_)?pathogenic/ && cv !~ /conflicting/)
    ben = (k ~ /[Bb]enign/)
    ours_p = (k ~ /athogenic/)
    if ((plp && (ben || k ~ /^VUS/)) || (ours_p && cv ~ /benign/)) {
        n++
        printf "  %-10s %-14s ours=%-22s ClinVar=%s (%s*)  %s\n", \
            $3, ($10 != "." ? $10 : substr($11, 1, 14)), $47 "(" $48 ")", substr($16, 1, 24), $17, $49
    }
}
END { if (!n) print "  none" }' "${results}/clinical_tiers.tsv"
echo ""
echo "10. LIMITATIONS"
echo "--------------"
echo "   - Short-read WGS: no CNV, structural variant, repeat-expansion or"
echo "     methylation analysis. Homologous / low-complexity regions unreliable."
echo "   - Mitochondrial variants are called with a diploid model (HaplotypeCaller);"
echo "     heteroplasmy is not quantified. Flagged chrM_diploid_caller."
echo "   - No phenotype (HPO) input: candidate tiers are not phenotype-prioritised."
echo "   - ACMG/AMP classification is automated and partial. Implemented codes:"
echo "     PVS1 (NMD-escape aware), PM2_Supporting, PM4, PP2, PP3/BP4 (ClinGen-"
echo "     calibrated REVEL and SpliceAI), BA1, BS1, BP7. Codes needing case-level"
echo "     or functional data are NOT evaluated: PS1-PS4, PM1, PM3, PM5, PM6, PP1,"
echo "     BS2-BS4, BP1-BP3, BP5. A VUS here therefore means insufficient automated"
echo "     evidence, not that the variant has been ruled out. Every classification"
echo "     requires review by a qualified variant scientist before clinical use."
echo "   - BA1 / BS1 use gnomAD v4.1 grpmax FAF95, retrieved per variant from the"
echo "     joint release. Rows whose Freq_Source still reads gnomADv2.1 could not"
echo "     be retrieved and their flags remain provisional. Freq_Max = NA means the"
echo "     variant was queried and is absent from gnomAD (PM2-supporting)."
echo "   - Gene-disease validity and inheritance come from GenCC (CC0). Genes with"
echo "     no GenCC assertion carry Gene_Validity = . and rank below curated genes;"
echo "     absence of an assertion is not evidence against the gene."
echo "   - Carrier status covers genes with a curated recessive or X-linked"
echo "     assertion; it is not an ACMG carrier-screening panel."
echo "   - No pharmacogenomic (PGx) interpretation."
echo "   - SpliceAI scores cover SNVs only (masked MANE v1.4 release)."
echo "   - QC_Flags / Confidence are computed from GATK annotations only; LOW"
echo "     confidence calls require orthogonal confirmation before any use."
echo ""
echo "11. DATA SOURCES"
echo "----------------"
echo "   Reference   : GRCh38 (hg38), MANE Select preferred transcripts"
echo "   Annotation  : GATK Funcotator (Gencode v43)"
echo "   ClinVar     : see funcotator data sources config"
echo "   Frequency   : gnomAD v4.1 joint (grpmax FAF95, remote lookup) + v2.1 fallback"
echo "   Predictors  : REVEL, AlphaMissense, SpliceAI (ClinGen-calibrated thresholds)"
echo "   Gene list   : ACMG SF v3.3 (84 genes)"
echo "   Gene-disease: GenCC submissions export (CC0 1.0)"
echo ""
echo "12. NOTE"
echo "--------"
echo "   Research / informational output. Not a diagnostic report. Any finding"
echo "   intended for clinical use requires orthogonal confirmation in an"
echo "   accredited laboratory and interpretation by a qualified clinician."
echo "==========================================================================="
} > "${REPORT}"

# ===========================================================================
# PART 4 - HTML report (self-contained, printable)
# ===========================================================================

echo "Writing HTML report..."
HTML="${results}/clinical_report.html"

{
cat << 'HTMLHEAD'
<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>GeneMolX Clinical Variant Report</title>
<style>
 :root{--bg:#fff;--fg:#1a1d23;--mut:#5b6472;--line:#e3e6eb;--accent:#0b5fff;--warn:#b45309;--bad:#b91c1c;--ok:#047857}
 *{box-sizing:border-box}
 body{margin:0;padding:0 16px 64px;font:14px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;color:var(--fg);background:var(--bg)}
 .wrap{max-width:1180px;margin:0 auto}
 header{padding:28px 0 18px;border-bottom:2px solid var(--fg)}
 h1{margin:0 0 4px;font-size:22px;letter-spacing:-.01em}
 .sub{color:var(--mut);font-size:13px}
 h2{margin:34px 0 10px;font-size:16px;border-bottom:1px solid var(--line);padding-bottom:6px}
 .cards{display:flex;flex-wrap:wrap;gap:10px;margin:18px 0}
 .card{flex:1 1 150px;border:1px solid var(--line);border-radius:8px;padding:10px 12px}
 .card .n{font-size:22px;font-weight:650}
 .card .l{font-size:12px;color:var(--mut)}
 table{width:100%;border-collapse:collapse;font-size:12.5px}
 th,td{text-align:left;padding:6px 8px;border-bottom:1px solid var(--line);vertical-align:top}
 th{font-weight:600;color:var(--mut);font-size:11.5px;text-transform:uppercase;letter-spacing:.04em}
 td.g{font-weight:600}
 .tag{display:inline-block;padding:1px 6px;border-radius:10px;font-size:11px;border:1px solid var(--line)}
 .HIGH{color:var(--ok);border-color:var(--ok)}
 .MEDIUM{color:var(--warn);border-color:var(--warn)}
 .LOW{color:var(--bad);border-color:var(--bad)}
 .flag{color:var(--bad);font-weight:600}
 .scroll{max-height:520px;overflow:auto;border:1px solid var(--line);border-radius:8px}
 .note{background:#f7f8fa;border-left:3px solid var(--accent);padding:10px 12px;margin:12px 0;font-size:13px}
 ul{margin:8px 0 0 18px;padding:0}li{margin:3px 0}
 footer{margin-top:40px;padding-top:14px;border-top:1px solid var(--line);color:var(--mut);font-size:12px}
 @media print{.scroll{max-height:none;overflow:visible}body{padding:0}}
</style></head><body><div class="wrap">
HTMLHEAD

if [ "${qc_pass}" -eq 0 ]; then printf '<div class="note" style="border-left-color:#b91c1c"><b>Sample failed quality control - these results are not interpretable.</b><br>%s<br>No variant has been classified. At this depth the variants passing the depth filter are pile-ups in collapsed repeats, which resemble loss-of-function calls but are artefacts.</div>\n' "${qc_reasons}"; fi
printf '<header><h1>Clinical Variant Report</h1><div class="sub">Sample <b>%s</b> &middot; generated %s &middot; GRCh38 &middot; GATK4 + Funcotator</div></header>\n' \
    "${sample}" "$(date '+%Y-%m-%d %H:%M')"

printf '<div class="cards">'
for t in "4:Secondary findings" "1:Known pathogenic" "5:Carrier" "2:Rare LoF" "3:Predicted damaging" "6:For review"; do
    printf '<div class="card"><div class="n">%s</div><div class="l">Tier %s &middot; %s</div></div>' \
        "$(cnt "${t%%:*}")" "${t%%:*}" "${t#*:}"
done
printf '</div>\n'

printf '<h2>Quality control</h2><table><tr><th>Metric</th><th>Value</th></tr>'
printf '<tr><td>Mean coverage</td><td>%s</td></tr>' "${mcov}"
printf '<tr><td>Genome &ge;10x / &ge;20x</td><td>%s / %s</td></tr>' "${p10}" "${p20}"
printf '<tr><td>Duplication rate</td><td>%s</td></tr>' "${dup}"
printf '<tr><td>Contamination (FREEMIX)</td><td>%s</td></tr>' "${freemix}"
printf '<tr><td>Reads aligned</td><td>%s</td></tr>' "${aln}"
printf '<tr><td>Genetic sex check</td><td>%s</td></tr>' "${sexline}"
printf '</table>\n'

html_rows() {
    mawk -F'\t' -v t="$1" '
    function esc(x) { gsub(/&/, "\\&amp;", x); gsub(/</, "\\&lt;", x); gsub(/>/, "\\&gt;", x); return x }
    BEGIN {
        print "<div class=\"scroll\"><table><tr><th>Gene</th><th>Variant</th><th>HGVS</th><th>Zyg</th>" \
              "<th>Effect</th><th>ACMG/AMP</th><th>Evidence codes</th><th>ClinVar</th><th>gnomAD v4.1</th>" \
              "<th>Inheritance</th><th>Gene validity</th><th>Conf</th><th>Quality flags</th><th>Why</th></tr>"
    }
    NR > 1 && $1 == t {
        n++
        f = ($24 == "NA") ? "not observed" : $24
        if ($26 != "." && $26 != "PM2_absent") f = f " <span class=\"flag\">" $26 "</span>"
        cl = $47; sub(/\(low_quality_call\)/, "", cl)
        printf "<tr><td class=\"g\">%s</td><td>%s:%s %s&gt;%s</td><td>%s</td><td>%s</td><td>%s</td>" \
               "<td><b>%s</b> (%s pts)</td><td>%s</td>" \
               "<td>%s (%s&#9733;)</td><td>%s</td><td>%s</td><td>%s</td>" \
               "<td><span class=\"tag %s\">%s</span></td><td>%s</td><td>%s</td></tr>\n", \
            esc($3), $4, $5, esc(substr($6,1,12)), esc(substr($7,1,12)), \
            esc(($10 != "." ? $10 : $11)), $8, esc($9), esc(cl), $48, esc($49), \
            esc($16), $17, f, esc($22), esc($23), \
            $30, $30, esc($31), esc($14)
    }
    END {
        if (!n) print "<tr><td colspan=\"14\">none</td></tr>"
        print "</table></div>"
    }' "${results}/clinical_tiers.tsv"
}

for t in 4 1 5 2 3 6; do
    case "$t" in
        4) h="Tier 4 &middot; Secondary findings (ACMG SF v3.3)";;
        1) h="Tier 1 &middot; Known pathogenic (ClinVar P/LP)";;
        5) h="Tier 5 &middot; Carrier status";;
        2) h="Tier 2 &middot; Rare loss-of-function candidates";;
        3) h="Tier 3 &middot; Rare predicted-damaging candidates";;
        6) h="Tier 6 &middot; Flagged for review";;
    esac
    printf '<h2>%s</h2>\n' "${h}"
    [ "$t" = "5" ] && [ "${plp_no_inh}" -gt 0 ] && printf '<div class="note">%s heterozygous P/LP variant(s) have no mode of inheritance in ACMG SF v3.3 or GenCC, so carrier status could not be assigned. It is not excluded for them.</div>\n' "${plp_no_inh}"
    html_rows "$t"
done

cat << 'HTMLTAIL'
<h2>Limitations</h2>
<ul>
 <li>Short-read WGS: no copy-number, structural variant, repeat-expansion or methylation analysis.</li>
 <li>Mitochondrial variants are called with a diploid model; heteroplasmy is not quantified.</li>
 <li>No phenotype (HPO) input, so candidate tiers are not phenotype-prioritised.</li>
 <li>Variants are not classified by ACMG/AMP criteria. Tier 1 reflects the ClinVar submitted classification.</li>
 <li>BA1 / BS1 use gnomAD v4.1 grpmax FAF95 where a record was retrieved; rows marked <i>_provisional</i> fall back to gnomAD v2.1 popmax.</li>
 <li>Carrier status covers genes with a curated inheritance assertion only; it is not an ACMG carrier-screening panel.</li>
 <li>No pharmacogenomic interpretation. SpliceAI covers SNVs only.</li>
 <li>Quality flags derive from GATK annotations; LOW-confidence calls require orthogonal confirmation.</li>
</ul>
<h2>Data sources</h2>
<ul>
 <li>Reference GRCh38 (hg38), MANE Select preferred transcripts</li>
 <li>GATK Funcotator, Gencode v43; ClinVar per the Funcotator data-source configuration</li>
 <li>gnomAD v4.1 joint release (remote FAF95 lookup) and gnomAD v2.1 (Funcotator)</li>
 <li>REVEL, AlphaMissense, SpliceAI at ClinGen-calibrated thresholds</li>
 <li>ACMG SF v3.3 (84 genes); GenCC gene-disease validity and inheritance (CC0 1.0)</li>
</ul>
<footer>
 Research / informational output. Not a diagnostic report. Any finding intended for clinical use
 requires orthogonal confirmation in an accredited laboratory and interpretation by a qualified clinician.
 <br>&copy; 2026 GeneMolX AI LLC. All rights reserved.
</footer>
</div></body></html>
HTMLTAIL
} > "${HTML}"

echo "  ${HTML}"

echo ""
echo "=== Clinical Relevant Variants Summary ==="
echo "Total P/LP (v1): $(( snp_count + indel_count )) variants (SNPs: ${snp_count}, INDELs: ${indel_count})"
echo "Tiered rows (v2): ${tier_total}"
echo "Output: ${results}/clinvar_plp_*.tsv, clinical_tiers.tsv, clinical_report.txt, clinical_report.html"
echo "Finished: $(date)"
