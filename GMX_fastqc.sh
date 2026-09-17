#!/bin/bash
# © 2026 GeneMolX AI LLC. All rights reserved.
# 1500 N GRANT ST STE N DENVER, CO, 80203, USA

set -euo pipefail

# Read quality control for the GeneMolX WGS pipeline - FIRST STEP
#
# Purpose: validate the raw data before anything else is done with it, and
# report the facts that decide the trimming parameters. Run it again on the
# trimmed reads to confirm the trimming worked.
#
# Accepts FASTQ as delivered by sequencing providers: a directory (searched
# recursively, so per-lane subfolders work), or individual files. Recognised
# extensions: .fastq.gz .fq.gz .fastq.bz2 .fq.bz2 .fastq .fq (any case).
# Recognised pair naming: _R1/_R2, _R1_001/_R2_001, _1/_2, .R1./.R2.
#
# --nogroup keeps every cycle separate in the per-base plots; FastQC otherwise
# bins cycles in reads over 50 bp and hides position-specific problems.

GMX_VERSION="0.1.0"
case "${1:-}" in
	-V|--version) echo "GMX_fastqc.sh ${GMX_VERSION} - GeneMolX WGS pipeline"; exit 0 ;;
esac

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate genemolx-wgs

echo "Started: $(date)"

usage() {
    cat << 'EOF'
Usage: GMX_fastqc.sh <output_dir> <reads_dir | fastq ...> [options]

  <output_dir>          where the reports are written
  <reads_dir | fastq>   a directory to search recursively, or one or more FASTQ files

Options:
  --threads N           FastQC threads (default 4)
  --include DIR         extra directory for the MultiQC report (e.g. aligned_reads);
                        may be given more than once
  --verify              full decompression integrity check of every input
                        (slow: minutes per file on a 30x genome)
  --preflight-reads N   records inspected for the format/encoding check (default 100000)
EOF
}

[ $# -ge 2 ] || { usage; exit 1; }

outdir="$1"; shift
threads=4
verify=0
pre_n=100000
inputs=()
include=()

while [ $# -gt 0 ]; do
    case "$1" in
        --threads)         threads="$2"; shift 2 ;;
        --include)         include+=("$2"); shift 2 ;;
        --verify)          verify=1; shift ;;
        --preflight-reads) pre_n="$2"; shift 2 ;;
        -h|--help)         usage; exit 0 ;;
        *)                 inputs+=("$1"); shift ;;
    esac
done

[ "${#inputs[@]}" -gt 0 ] || { echo "ERROR: no input given"; usage; exit 1; }
mkdir -p "${outdir}"

# -------------------
# STEP 1: Resolve inputs
# -------------------

fastqs=()
for item in "${inputs[@]}"; do
    if [ -d "${item}" ]; then
        while IFS= read -r f; do fastqs+=("${f}"); done < <(find "${item}" -type f \
            \( -iname '*.fastq.gz' -o -iname '*.fq.gz' -o -iname '*.fastq.bz2' \
               -o -iname '*.fq.bz2' -o -iname '*.fastq' -o -iname '*.fq' \) | sort)
        # a provider that ships .zip cannot be read by FastQC directly
        if find "${item}" -type f -iname '*.fastq.zip' -o -iname '*.fq.zip' | command grep -q .; then
            echo "NOTE: .zip FASTQ found - FastQC cannot read zip; decompress or re-gzip first"
        fi
    elif [ -f "${item}" ]; then
        fastqs+=("${item}")
    else
        echo "ERROR: not found: ${item}"; exit 1
    fi
done

[ "${#fastqs[@]}" -gt 0 ] || { echo "ERROR: no FASTQ files found"; exit 1; }

# -------------------
# STEP 2: Preflight - is this really usable FASTQ, and how is it encoded?
# -------------------

echo ""
echo "=== Input validation (${#fastqs[@]} file(s)) ==="
printf "  %-42s %-10s %-9s %-11s %s\n" FILE SIZE ENCODING READ_LEN PLATFORM

PRE_AWK='
NR % 4 == 1 { if (substr($0,1,1) != "@") { print "BADHDR"; exit } if (NR == 1) hdr = $0 }
NR % 4 == 2 { len = length($0); if (len < minl || minl == 0) minl = len; if (len > maxl) maxl = len; sl = len }
NR % 4 == 3 { if (substr($0,1,1) != "+") { print "BADSEP"; exit } }
NR % 4 == 0 {
    if (length($0) != sl) { print "LENMISMATCH"; exit }
    for (i = 1; i <= length($0); i++) { c = index(Q, substr($0,i,1)) + 31; if (c < lo || lo == 0) lo = c; if (c > hi) hi = c }
    n++
    if (n >= maxn) exit
}
END {
    if (n == 0) { print "EMPTY"; exit }
    enc = (lo < 59) ? "phred33" : ((hi > 74) ? "phred64" : "phred33?")
    # platform guess from the read-name convention
    nc = split(hdr, h, ":")
    if (nc >= 7)                         plat = "Illumina(" substr(h[1], 2) ")"
    else if (hdr ~ /\/[12]$/ && hdr ~ /L[0-9]+C[0-9]+R[0-9]+/) plat = "MGI/DNBSEQ"
    else if (nc >= 5)                    plat = "Illumina-like"
    else                                 plat = "unknown"
    print "OK", enc, minl "-" maxl, plat, lo, hi, n
}'
Q=' !"#$%&'"'"'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'

decomp() {  # stream a file whatever its compression
    case "$1" in
        *.gz|*.GZ)   zcat "$1" ;;
        *.bz2|*.BZ2) bzcat "$1" ;;
        *)           cat "$1" ;;
    esac
}

fail=0
declare -A ENC
for f in "${fastqs[@]}"; do
    size=$(du -h "$f" | cut -f1)
    res=$(set +o pipefail; decomp "$f" 2> /dev/null | head -n $(( pre_n * 4 )) \
          | mawk -v Q="${Q}" -v maxn="${pre_n}" "${PRE_AWK}")
    set -- ${res}
    case "${1:-EMPTY}" in
        OK) printf "  %-42s %-10s %-9s %-11s %s\n" "$(basename "$f")" "${size}" "$2" "$3" "$4"
            ENC["$f"]="$2" ;;
        BADHDR)      echo "  FAIL $(basename "$f"): record header does not start with @"; fail=1 ;;
        BADSEP)      echo "  FAIL $(basename "$f"): third line of a record does not start with +"; fail=1 ;;
        LENMISMATCH) echo "  FAIL $(basename "$f"): sequence and quality lengths differ"; fail=1 ;;
        EMPTY)       echo "  FAIL $(basename "$f"): no readable FASTQ records (corrupt or wrong format)"; fail=1 ;;
        *)           echo "  FAIL $(basename "$f"): unreadable"; fail=1 ;;
    esac

    if [ "${verify}" -eq 1 ]; then
        case "$f" in
            *.gz|*.GZ)   gzip -t "$f"  2> /dev/null || { echo "  FAIL $(basename "$f"): gzip integrity"; fail=1; } ;;
            *.bz2|*.BZ2) bzip2 -t "$f" 2> /dev/null || { echo "  FAIL $(basename "$f"): bzip2 integrity"; fail=1; } ;;
        esac
    fi
done
[ "${fail}" -eq 0 ] || { echo ""; echo "ERROR: input validation failed - fix the files above before continuing"; exit 1; }

# -------------------
# STEP 3: Pair detection
# -------------------

echo ""
echo "=== Pairing ==="
for f in "${fastqs[@]}"; do basename "$f"; done | mawk '
function key(n,   k) {
    # replace only the mate marker, keeping any suffix, so _R1_paired pairs with
    # _R2_paired and not with _R1_unpaired
    k = n; if (sub(/_R1([._])/, "_@@\\1", k)) return k "\t1"
    k = n; if (sub(/_R2([._])/, "_@@\\1", k)) return k "\t2"
    k = n; if (sub(/\.R1\./, ".@@.", k))      return k "\t1"
    k = n; if (sub(/\.R2\./, ".@@.", k))      return k "\t2"
    k = n; if (sub(/_1\.f(ast)?q/, "_@@.fq", k)) return k "\t1"
    k = n; if (sub(/_2\.f(ast)?q/, "_@@.fq", k)) return k "\t2"
    return n "\t0"
}
{ split(key($0), a, "\t"); if (a[2] == "1") r1[a[1]] = $0; else if (a[2] == "2") r2[a[1]] = $0; else single[$0] = 1 }
END {
    for (b in r1) if (b in r2) { printf "  pair    %s  +  %s\n", r1[b], r2[b]; delete r2[b] }
                  else           printf "  R1 only %s  (mate not found)\n", r1[b]
    for (b in r2) printf "  R2 only %s  (mate not found)\n", r2[b]
    for (s in single) printf "  single  %s\n", s
}'

# -------------------
# STEP 4: FastQC
# -------------------

echo ""
echo "FastQC: ${#fastqs[@]} file(s), ${threads} thread(s)"
fastqc --nogroup --memory 2048 -t "${threads}" -o "${outdir}" "${fastqs[@]}"

# -------------------
# STEP 5: Aggregate report
# -------------------

if command -v multiqc > /dev/null 2>&1; then
    echo "MultiQC: aggregating"
    multiqc --quiet --force --outdir "${outdir}" "${outdir}" ${include[@]+"${include[@]}"}
    echo "  ${outdir}/multiqc_report.html"
else
    echo "MultiQC not installed - per-file FastQC reports only"
    echo "  install it with: conda install -n genemolx-wgs -c bioconda multiqc"
fi

# -------------------
# STEP 6: Verdicts and trimming guidance
# -------------------

echo ""
echo "=== FastQC module verdicts (WARN / FAIL only) ==="
found=0
for z in "${outdir}"/*_fastqc.zip; do
    [ -e "${z}" ] || continue
    name="$(basename "${z}" _fastqc.zip)"
    if unzip -p "${z}" "*/summary.txt" 2> /dev/null | mawk -F'\t' -v n="${name}" \
        '$1 != "PASS" { printf "  %-8s %-40s %s\n", $1, $2, n; f = 1 } END { exit !f }'; then
        found=1
    fi
done
[ "${found}" -eq 1 ] || echo "  none - all modules passed"

echo ""
echo "=== Guidance for GMX_trim_reads.sh ==="
{
    for z in "${outdir}"/*_fastqc.zip; do
        [ -e "${z}" ] || continue
        unzip -p "${z}" "*/fastqc_data.txt" 2> /dev/null
    done
} | mawk -F'\t' '
function worst(cur, new) { if (new == "fail" || cur == "fail") return "fail"
                           if (new == "warn" || cur == "warn") return "warn"
                           return (new == "" ? cur : new) }
/^Encoding\t/                      { enc[$2] = 1 }
/^Sequence length\t/               { len[$2] = 1
                                     l = $2; sub(/^[0-9]+-/, "", l)
                                     if (l + 0 > maxlen) maxlen = l + 0 }
/^>>Adapter Content\t/             { ad = worst(ad, $2) }
/^>>Overrepresented sequences\t/   { ov = worst(ov, $2) }
/^Total Sequences\t/               { tot += $2 }
END {
    printf "  reads (all files)      : %d\n", tot
    es = ""; for (e in enc) es = es e "; "; printf "  quality encoding       : %s\n", es
    ls = ""; for (l in len) ls = ls l " ";  printf "  read length            : %s\n", ls
    printf "  adapter content        : %s\n", (ad == "" ? "n/a" : ad)
    printf "  overrepresented seqs   : %s\n", (ov == "" ? "n/a" : ov)
    print  ""
    if (es ~ /Illumina 1\.[35]/ || es ~ /Solexa/)
        print "  ! phred64 encoding detected - pass -phred64 to Trimmomatic"
    if (ad == "fail" || ad == "warn")
        print "  ! adapter contamination present - keep ILLUMINACLIP (resources/adapters/combined_adapters.fa)"
    else
        print "  adapter content is clean; ILLUMINACLIP can stay as a safeguard"
    if (maxlen > 0) {
        m = int(maxlen / 3)
        m = int((m + 2) / 5) * 5          # round to the nearest 5
        if (m < 36) m = 36                 # below ~36 bp a read maps ambiguously
        printf "  suggested MINLEN       : %d   (reads are %d bp; MINLEN ~ 1/3 of read length)\n", m, maxlen
        if (m != 50) printf "  ! the script default is MINLEN:50 - pass %d for this data\n", m
    } else
        print  "  set MINLEN below the shortest read you want to keep (default 50 suits 150 bp reads)"
}'

echo ""
echo "Reports: ${outdir}"
echo "Finished: $(date)"
