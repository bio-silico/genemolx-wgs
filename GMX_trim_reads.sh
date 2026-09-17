#!/bin/bash
# © 2026 GeneMolX AI LLC. All rights reserved.
# 1500 N GRANT ST STE N DENVER, CO, 80203, USA

set -euo pipefail

# Activate conda environment
GMX_VERSION="0.1.0"
case "${1:-}" in
	-V|--version) echo "GMX_trim_reads.sh ${GMX_VERSION} - GeneMolX WGS pipeline"; exit 0 ;;
esac

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate genemolx-wgs

# Usage: GMX_trim_reads.sh <R1.fq.gz> <R2.fq.gz> <output_dir> <sample_id> [adapter_file] [threads] [phred]
#   phred: auto (default) | phred33 | phred64
#   auto detects the quality encoding from the first reads of R1 and passes it to
#   Trimmomatic explicitly, so the encoding used is recorded in the log instead of
#   being inferred silently.
if [ $# -lt 4 ]; then
	echo "Usage: GMX_trim_reads.sh <R1.fq.gz> <R2.fq.gz> <output_dir> <sample_id> [adapter_file] [threads] [phred]"
	echo "  phred: auto (default) | phred33 | phred64"
	exit 1
fi

R1="$1"
R2="$2"
output_dir="$3"
sample="$4"
adapter_file="${5:-$(dirname "$0")/resources/adapters/combined_adapters.fa}"
threads="${6:-4}"
phred="${7:-auto}"

for f in "${R1}" "${R2}"; do
	[ -f "${f}" ] || { echo "ERROR: input FASTQ not found: ${f}"; exit 1; }
done
[ -f "${adapter_file}" ] || { echo "ERROR: adapter file not found: ${adapter_file}"; exit 1; }

mkdir -p "${output_dir}"

# Detect the quality encoding from the first 100k reads of R1 (same test as GMX_fastqc.sh)
if [ "${phred}" = "auto" ]; then
	Q=' !"#$%&'"'"'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'
	case "${R1}" in
		*.gz|*.GZ)   reader=zcat ;;
		*.bz2|*.BZ2) reader=bzcat ;;
		*)           reader=cat ;;
	esac
	phred=$(set +o pipefail; ${reader} "${R1}" 2> /dev/null | head -n 400000 | mawk -v Q="${Q}" '
		NR % 4 == 0 {
			for (i = 1; i <= length($0); i++) { c = index(Q, substr($0, i, 1)) + 31
				if (c < lo || lo == 0) lo = c; if (c > hi) hi = c }
			n++
		}
		END { print (n == 0) ? "phred33" : ((lo < 59) ? "phred33" : ((hi > 74) ? "phred64" : "phred33")) }')
	echo "Quality encoding detected: ${phred}"
fi

# Define output files
R1_paired="${output_dir}/${sample}_R1_paired.fastq.gz"
R1_unpaired="${output_dir}/${sample}_R1_unpaired.fastq.gz"
R2_paired="${output_dir}/${sample}_R2_paired.fastq.gz"
R2_unpaired="${output_dir}/${sample}_R2_unpaired.fastq.gz"

# Run Trimmomatic with updated adapter sequences and additional trimming options
trimmomatic PE -threads "${threads}" -"${phred}" \
    "${R1}" "${R2}" \
    "${R1_paired}" "${R1_unpaired}" \
    "${R2_paired}" "${R2_unpaired}" \
    ILLUMINACLIP:"${adapter_file}":2:30:10:8:true \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 AVGQUAL:20 \
    2> "${output_dir}/${sample}_trimmomatic.log"

# Print a completion message
echo "Trimming completed successfully!"
