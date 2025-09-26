#!/bin/bash

# --- Stricter Shell Safety ---
set -e
set -o pipefail
IFS=$'\n\t'

# Description:
#   [cite_start]This script executes a csRNA-seq workflow based on the HOMER tutorial[cite: 1, 2].
#   It processes csRNA, input, and RNA-seq (SE or PE) samples, performs
#   TSS calling, differential expression, and motif analysis.

# --- 1. Default Parameters & Help Function ---
THREADS=8
SIF_PATH="$(dirname "$0")/csRNA.sif"

HELP_MSG="Usage: $0 -s <sample_sheet.csv> -o <out_dir> -r <ref_dir> -c <container.sif> [OPTIONS]

Required:
  -s  Sample sheet (CSV). Required 6-column format:
      sample,condition,type,fastq_r1,fastq_r2,strandedness
      (leave fastq_r2 blank for single-end samples)
  -o  Main output directory.
  -r  Reference data directory (must contain STAR_index/, rRNA_STAR_index/, a .fa, and a .gtf file).
  -c  Path to the Singularity container (csRNA.sif).

Optional:
  -t  Number of threads to use (Default: 8).
  -h  Display this help message.
"

# --- 2. Parse Command-line Arguments ---
while getopts "s:o:r:c:t:h" opt; do
  case ${opt} in
    s ) SAMPLE_SHEET=$(realpath "${OPTARG}") ;;
    o ) OUT_DIR=$(realpath "${OPTARG}") ;;
    r ) REF_DIR=$(realpath "${OPTARG}") ;;
    c ) SIF_PATH=$(realpath "${OPTARG}") ;;
    t ) THREADS=${OPTARG} ;;
    h ) echo "${HELP_MSG}"; exit 0 ;;
    \? ) echo "Invalid option: -${OPTARG}" >&2; echo "${HELP_MSG}"; exit 1 ;;
  esac
done

if [ -z "${SAMPLE_SHEET}" ] || [ -z "${OUT_DIR}" ] || [ -z "${REF_DIR}" ] || [ -z "${SIF_PATH}" ]; then
    echo "Error: Missing mandatory arguments." >&2; echo "${HELP_MSG}"; exit 1
fi

# --- 3. Setup Environment and Directories ---
mkdir -p "${OUT_DIR}"/{qc/{raw_fastqc,trimmed_fastqc,trimmed_reads},alignments,alignments/rRNA_logs,tag_directories,visualization,tss_calling,differential_expression,motif_analysis}
GENOME_FASTA=$(find "${REF_DIR}" -name "GRCh38.primary_assembly.genome.fa" | head -n 1)
GTF_FILE=$(find "${REF_DIR}" -name "annotation.gtf" | head -n 1)
GENOME_INDEX_DIR="${REF_DIR}/star_genome_index"
RRNA_INDEX_DIR="${REF_DIR}/star_rrna_index" # Path for rRNA STAR index

if [ -z "${GENOME_FASTA}" ] || [ -z "${GTF_FILE}" ]; then
    echo "Error: Could not find genome .fa and .gtf files in ${REF_DIR}" >&2; exit 1
fi
if [ ! -d "${GENOME_INDEX_DIR}" ]; then
    echo "Error: STAR genome index not found at ${GENOME_INDEX_DIR}" >&2; exit 1
fi

echo "================================================="
echo "====== csRNA-seq HOMER Pipeline Started ======="
echo "================================================="
echo "Sample Sheet: ${SAMPLE_SHEET}"
echo "Output Directory: ${OUT_DIR}"
echo "Reference Directory: ${REF_DIR}"
echo "Threads: ${THREADS}"
echo "================================================="

# --- 4. Parse Sample Sheet and Process Each File ---
declare -A SAMPLES
# Read sample sheet (6 columns), skipping header, and populate the associative array
while IFS=, read -r sample condition type fastq_r1 fastq_r2 strandedness; do
    if [ -z "$sample" ] || [ -z "$type" ]; then continue; fi
    SAMPLES["${sample}_${type}_path"]="$fastq_r1"
    SAMPLES["${sample}_condition"]="$condition"
    if [ "$type" == "rna" ]; then
        SAMPLES["${sample}_rna_path2"]="$fastq_r2"
        SAMPLES["${sample}_rna_strandedness"]="$strandedness"
    fi
done < <(tail -n +2 "${SAMPLE_SHEET}" | tr -d '\r')

UNIQUE_SAMPLES=$(tail -n +2 "${SAMPLE_SHEET}" | tr -d '\r' | cut -d, -f1 | sort -u)

for sample in ${UNIQUE_SAMPLES}; do
    echo -e "\n--- Processing sample group: [${sample}] ---"
    
    # --- Process csRNA and input files ---
    for type in csRNA input; do
        fastq_path_key="${sample}_${type}_path"
        if [[ -v SAMPLES[$fastq_path_key] ]] && [ -n "${SAMPLES[$fastq_path_key]}" ]; then
            FQ_PATH=${SAMPLES[$fastq_path_key]}
            FQ_NAME=$(basename "${FQ_PATH}" .fastq.gz)
            echo "--- Found ${type} file: ${FQ_PATH}"

            TRIMMED_FQ_OUTDIR="${OUT_DIR}/qc/trimmed_reads"
            TRIMMED_FQ="${TRIMMED_FQ_OUTDIR}/${FQ_NAME}_trimmed.fq.gz"
            ALIGNED_SAM="${OUT_DIR}/alignments/${FQ_NAME}.Aligned.out.sam"
            TAG_DIR="${OUT_DIR}/tag_directories/${sample}_${type}"

            echo "[${sample}/${type}] Running FastQC -> Trim Galore -> FastQC..."
            singularity exec "${SIF_PATH}" fastqc "${FQ_PATH}" -t "${THREADS}" -o "${OUT_DIR}/qc/raw_fastqc"
            singularity exec "${SIF_PATH}" trim_galore --cores "${THREADS}" --gzip -o "${TRIMMED_FQ_OUTDIR}" --stringency 4 --length 20 "${FQ_PATH}"
            singularity exec "${SIF_PATH}" fastqc "${TRIMMED_FQ}" -t "${THREADS}" -o "${OUT_DIR}/qc/trimmed_fastqc"

            echo "[${sample}/${type}] Aligning with STAR..."
            singularity exec -B "${GENOME_INDEX_DIR}:/genome" "${SIF_PATH}" STAR \
                --genomeDir /genome --runThreadN "${THREADS}" --readFilesIn "${TRIMMED_FQ}" --readFilesCommand zcat \
                --outFileNamePrefix "${OUT_DIR}/alignments/${FQ_NAME}." --outSAMstrandField intronMotif \
                --outMultimapperOrder Random --outSAMmultNmax 1 --outFilterMultimapNmax 10000 --limitOutSAMoneReadBytes 10000000
            
            echo "[${sample}/${type}] Creating Tag Directory..."
            singularity exec -B "${REF_DIR}:${REF_DIR}" "${SIF_PATH}" makeTagDirectory "${TAG_DIR}" "${ALIGNED_SAM}" -genome "${GENOME_FASTA}" -checkGC -fragLength 150
        fi
    done

    # --- Process RNA-seq file(s) ---
    type="rna"
    fastq_path_key="${sample}_${type}_path"
    if [[ -v SAMPLES[$fastq_path_key] ]] && [ -n "${SAMPLES[$fastq_path_key]}" ]; then
        FQ_PATH_R1=${SAMPLES[$fastq_path_key]}
        FQ_NAME_R1=$(basename "${FQ_PATH_R1}" .fastq.gz)
        
        # Check for paired-end data
        fastq_path_key_r2="${sample}_rna_path2"
        IS_PAIRED_END=false
        if [[ -v SAMPLES[$fastq_path_key_r2] ]] && [ -n "${SAMPLES[$fastq_path_key_r2]}" ]; then
            IS_PAIRED_END=true
            FQ_PATH_R2=${SAMPLES[$fastq_path_key_r2]}
            FQ_NAME_R2=$(basename "${FQ_PATH_R2}" .fastq.gz)
            echo "--- Found RNA-seq Paired-End data for ${sample}"
        else
            echo "--- Found RNA-seq Single-End data for ${sample}"
        fi

        TRIMMED_FQ_OUTDIR="${OUT_DIR}/qc/trimmed_reads"
        ALIGNED_SAM="${OUT_DIR}/alignments/${FQ_NAME_R1}.Aligned.out.sam"
        TAG_DIR="${OUT_DIR}/tag_directories/${sample}_${type}"

        # Step 1: QC and Trimming
        echo "[${sample}/rna] Running FastQC -> Trim Galore -> FastQC..."
        singularity exec "${SIF_PATH}" fastqc "${FQ_PATH_R1}" -t "${THREADS}" -o "${OUT_DIR}/qc/raw_fastqc"
        if [ "$IS_PAIRED_END" = true ]; then
            singularity exec "${SIF_PATH}" fastqc "${FQ_PATH_R2}" -t "${THREADS}" -o "${OUT_DIR}/qc/raw_fastqc"
            singularity exec "${SIF_PATH}" trim_galore --cores "${THREADS}" --gzip -o "${TRIMMED_FQ_OUTDIR}" --paired --length 20 "${FQ_PATH_R1}" "${FQ_PATH_R2}"
            TRIMMED_FQ_R1="${TRIMMED_FQ_OUTDIR}/${FQ_NAME_R1}_val_1.fq.gz"
            TRIMMED_FQ_R2="${TRIMMED_FQ_OUTDIR}/${FQ_NAME_R2}_val_2.fq.gz"
            singularity exec "${SIF_PATH}" fastqc "${TRIMMED_FQ_R1}" "${TRIMMED_FQ_R2}" -t "${THREADS}" -o "${OUT_DIR}/qc/trimmed_fastqc"
        else
            singularity exec "${SIF_PATH}" trim_galore --cores "${THREADS}" --gzip -o "${TRIMMED_FQ_OUTDIR}" --length 20 "${FQ_PATH_R1}"
            TRIMMED_FQ_R1="${TRIMMED_FQ_OUTDIR}/${FQ_NAME_R1}_trimmed.fq.gz"
            singularity exec "${SIF_PATH}" fastqc "${TRIMMED_FQ_R1}" -t "${THREADS}" -o "${OUT_DIR}/qc/trimmed_fastqc"
        fi

        # Step 2: Alignment (rRNA removal -> Genome alignment)
        echo "[${sample}/rna] Removing rRNA reads..."
        if [ ! -d "${RRNA_INDEX_DIR}" ]; then
            echo "Error: RNA-seq sample found, but rRNA STAR index not found at ${RRNA_INDEX_DIR}" >&2; exit 1
        fi
        RRNA_LOG_PREFIX="${OUT_DIR}/alignments/rRNA_logs/${FQ_NAME_R1}.rRNA."
        READS_IN_RRNA="${TRIMMED_FQ_R1}"
        [ "$IS_PAIRED_END" = true ] && READS_IN_RRNA+=" ${TRIMMED_FQ_R2}"
        
        singularity exec -B "${RRNA_INDEX_DIR}:/rRNA_genome" "${SIF_PATH}" STAR \
            --genomeDir /rRNA_genome --runThreadN "${THREADS}" --readFilesIn ${READS_IN_RRNA} --readFilesCommand zcat \
            --outFileNamePrefix "${RRNA_LOG_PREFIX}" --outReadsUnmapped Fastx
        
        UNMAPPED_R1="${RRNA_LOG_PREFIX}Unmapped.out.mate1"
        READS_IN_GENOME="${UNMAPPED_R1}"
        [ "$IS_PAIRED_END" = true ] && READS_IN_GENOME+=" ${RRNA_LOG_PREFIX}Unmapped.out.mate2"

        echo "[${sample}/rna] Aligning to genome with standard RNA-seq parameters..."
        singularity exec -B "${GENOME_INDEX_DIR}:/genome" "${SIF_PATH}" STAR \
            --genomeDir /genome --runThreadN "${THREADS}" --readFilesIn ${READS_IN_GENOME} \
            --outFileNamePrefix "${OUT_DIR}/alignments/${FQ_NAME_R1}." \
            --twopassMode Basic --outSAMstrandField intronMotif

        # Step 3: Create Tag Directory
        echo "[${sample}/rna] Creating Tag Directory..."
        strandedness_key="${sample}_rna_strandedness"
        strand_opt=${SAMPLES[$strandedness_key]}
        MAKEDIR_OPTS=""
        case ${strand_opt} in
            SE_plus)    MAKEDIR_OPTS="" ;;
            PE_plus)    MAKEDIR_OPTS="-read1" ;;
            SE_minus)   MAKEDIR_OPTS="-flip" ;;
            PE_minus)   MAKEDIR_OPTS="-read2" ;;
            *)          echo "Warning: Invalid strandedness '${strand_opt}' for ${sample}" >&2 ;;
        esac
        echo "Using RNA-seq strandedness option: ${MAKEDIR_OPTS}"
        singularity exec "${SIF_PATH}" makeTagDirectory "${TAG_DIR}" "${ALIGNED_SAM}" ${MAKEDIR_OPTS}
    fi

    # --- Step 4: Create Visualization Files (for csRNA and input) ---
    echo "[${sample}] Creating visualization files..."
    for type in csRNA input; do
        TAG_DIR_TO_VISUALIZE="${OUT_DIR}/tag_directories/${sample}_${type}"
        if [ -d "${TAG_DIR_TO_VISUALIZE}" ]; then
            echo "... for ${type}"
            singularity exec "${SIF_PATH}" makeUCSCfile "${TAG_DIR_TO_VISUALIZE}" -style tss -strand + > "${OUT_DIR}/visualization/${sample}_${type}.pos.bedGraph"
            singularity exec "${SIF_PATH}" makeUCSCfile "${TAG_DIR_TO_VISUALIZE}" -style tss -strand - -neg > "${OUT_DIR}/visualization/${sample}_${type}.neg.bedGraph"
        fi
    done
done

# --- 5. Finding TSS Clusters (Full Approach) ---
echo -e "\n\n--- Step 5: Finding TSS clusters for each sample ---"
for sample in ${UNIQUE_SAMPLES}; do
    CSRNA_TAG_DIR="${OUT_DIR}/tag_directories/${sample}_csRNA"
    INPUT_TAG_DIR="${OUT_DIR}/tag_directories/${sample}_input"
    RNA_TAG_DIR="${OUT_DIR}/tag_directories/${sample}_rna"
    
    if [ ! -d "${CSRNA_TAG_DIR}" ]; then
        echo "Skipping TSS calling for ${sample}, csRNA data not found."
        continue
    fi

    FINDTSS_CMD="findcsRNATSS.pl ${CSRNA_TAG_DIR} -o ${OUT_DIR}/tss_calling/${sample} -gtf ${GTF_FILE} -genome ${GENOME_FASTA}"
    [ -d "${INPUT_TAG_DIR}" ] && FINDTSS_CMD+=" -i ${INPUT_TAG_DIR}" || echo "Warning: No input control found for sample ${sample}."
    [ -d "${RNA_TAG_DIR}" ] && FINDTSS_CMD+=" -rna ${RNA_TAG_DIR}"

    echo "Running for ${sample}: ${FINDTSS_CMD}"
    eval "singularity exec -B "${REF_DIR}:${REF_DIR}" "${SIF_PATH}" ${FINDTSS_CMD}"
done

# --- 6. Differential Expression Analysis ---
echo -e "\n\n--- Step 6: Quantifying and Calculating Differential Expression ---"

MERGED_PEAKS="${OUT_DIR}/differential_expression/merged_tss_clusters.txt"
RAW_COUNTS="${OUT_DIR}/differential_expression/raw_counts.txt"
DIFF_RESULTS="${OUT_DIR}/differential_expression/differential_results.txt"

# Step 6.1 Merge all TSS clusters
echo "Merging TSS clusters from all samples..."
TSS_FILES=$(find "${OUT_DIR}/tss_calling" -name "*.tss.txt")
if [ -z "${TSS_FILES}" ]; then
    echo "No TSS files found to merge. Skipping downstream analysis."
    exit 0
fi
singularity exec "${SIF_PATH}" mergePeaks -strand ${TSS_FILES} > "${MERGED_PEAKS}"

# Step 6.2 Quantify csRNA read counts across merged clusters
echo "Quantifying reads in merged TSS clusters..."
for s in ${UNIQUE_SAMPLES}; do
    echo "${OUT_DIR}/tag_directories/${s}_csRNA/" >> tagDirs.txt
done

singularity exec "${SIF_PATH}" annotatePeaks.pl "${MERGED_PEAKS}" "${GENOME_FASTA}" \
    -strand + -fragLength 1 -raw -dfile tagDirs.txt > "${RAW_COUNTS}"

# Step 6.3 Run differential expression if replicates exist
if [ $(wc -l < "${SAMPLE_SHEET}") -gt 3 ]; then
    echo "Performing differential expression with DESeq2..."

singularity exec "${SIF_PATH}" Rscript - <<'EOF'
suppressMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
})

raw_counts_file <- Sys.getenv("RAW_COUNTS")
sample_sheet_file <- Sys.getenv("SAMPLE_SHEET")
output_file <- Sys.getenv("DIFF_RESULTS")

# 1. 读入 sample sheet
samples <- read.csv(sample_sheet_file, stringsAsFactors = FALSE)
samples <- samples[samples$type == "csRNA", c("sample","condition")]
rownames(samples) <- samples$sample

# 2. 读入 raw_counts
raw <- read.delim(raw_counts_file, header = TRUE, check.names = FALSE)
count_cols <- grep("Tag Count", colnames(raw), value = TRUE)

# 3. 构建 count 矩阵
counts <- raw[, count_cols]
rownames(counts) <- raw$PeakID
colnames(counts) <- sub(".*/([^/]+)/ Tag.*", "\\1", count_cols)

# 4. 确保样本顺序一致
counts <- counts[, rownames(samples)]

# 5. 构建 DESeq2 对象
dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(counts)),
                              colData = samples,
                              design = ~ condition)
dds <- DESeq(dds)

# 6. 自动检测条件组
conds <- levels(dds$condition)
if(length(conds) != 2){
  stop("Condition groups must be exactly 2 for DESeq2 contrast.")
}
res <- results(dds, contrast=c("condition", conds[2], conds[1]))
res <- as.data.frame(res)

# 7. 合并到原始注释信息
final <- cbind(raw[, 1:(ncol(raw)-ncol(counts))], counts, res)

# 8. 输出
write.table(final, file=output_file, sep="\t", quote=FALSE, row.names=FALSE)
EOF

else
    echo "Warning: Fewer than 2 replicates per condition – skipping differential expression."
fi


# --- 7. Motif Analysis ---
echo -e "\n\n--- Step 7: Motif Discovery and Sequence Feature Analysis ---"
MOTIF_DIR="${OUT_DIR}/motif_analysis"
mkdir -p "${MOTIF_DIR}"

# Step 7.1 Nucleotide frequency analysis around TSS
echo "Analyzing nucleotide composition near TSS..."
singularity exec "${SIF_PATH}" annotatePeaks.pl "${MERGED_PEAKS}" "${GENOME_FASTA}" \
    -size 1000 -hist 1 -di > "${MOTIF_DIR}/nucleotide_freq.txt"

# Step 7.2 Motif enrichment analysis
echo "Performing motif discovery..."
singularity exec "${SIF_PATH}" findMotifsGenome.pl "${MERGED_PEAKS}" "${GENOME_FASTA}" \
    "${MOTIF_DIR}" -size -150,50 -len 8,10,12 -p "${THREADS}"

# Step 7.3 Motif distribution relative to TSS
echo "Checking distribution of example motif (e.g., Sp1)..."
singularity exec "${SIF_PATH}" annotatePeaks.pl "${MERGED_PEAKS}" "${GENOME_FASTA}" \
    -size 500 -hist 1 -m "${MOTIF_DIR}/knownResults/known1.motif" > "${MOTIF_DIR}/motif_distribution.txt"

# --- 8. Final MultiQC Report ---
echo -e "\n\n--- Generating Final MultiQC Report ---"
eval "singularity exec "${SIF_PATH}" multiqc "${OUT_DIR}" -o "${OUT_DIR}/multiqc_report" --force"

echo -e "\n\n==============================================="
echo "====== csRNA-seq Pipeline Completed Successfully! ======"
echo "==============================================="
echo "Results are in: ${OUT_DIR}"
echo "View differential expression at: ${DIFF_RESULTS}"
echo "View motif analysis at: ${MOTIF_DIR}/homerResults.html"
echo "==============================================="