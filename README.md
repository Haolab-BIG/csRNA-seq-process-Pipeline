# csRNA-seq Pipeline

This unified csRNA-seq pipeline processes raw FASTQ files through to differential TSS (Transcription Start Site) analysis and motif discovery. It uses Singularity for reproducibility, incorporating csRNA-seq, input control, and optional total RNA-seq data (both single-end and paired-end).

## Workflow

<img width="2174" height="454" alt="CleanShot 2025-09-25 at 18 17 53@2x" src="https://github.com/user-attachments/assets/df8e60b7-f06b-4da4-b71a-4b78cd7071a5" />


## Features

  * **Single Command Execution**: Executes the entire workflow—from FASTQ QC, alignment, and TSS calling, through read quantification, to differential expression analysis and motif discovery—with a single shell command.
  * **Reproducible**: All core software (FastQC, Trim Galore, STAR, HOMER, R/DESeq2) is encapsulated within a Singularity container (csRNA.sif), ensuring the analysis is fully reproducible across different computing environments.
  * **Automated Reporting**: Generates a final, interactive MultiQC report summarizing quality control metrics across all samples and steps for easy assessment.

## Requirements

1.  **System Configuration**:
    * 8-core CPU
    * 64 GB RAM

2.  **Singularity**: Must be installed on your system. Below are detailed steps for installing on an Ubuntu 22.04 system. For other operating systems, please refer to the [official installation guide](https://www.google.com/search?q=https://docs.sylabs.io/guides/latest/user-guide/installation.html).

      * **Step 1: Install System Dependencies**

        ```bash
        # Update package lists and install dependencies
        sudo apt-get update
        sudo apt-get install -y \
            build-essential \
            libseccomp-dev \
        	libfuse3-dev \
            pkg-config \
            squashfs-tools \
            cryptsetup \
            curl wget git
        ```

      * **Step 2: Install Go Language**

        ```bash
        # Download and install Go (check for the latest version)
        wget https://go.dev/dl/go1.21.3.linux-amd64.tar.gz
        sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
        rm go1.21.3.linux-amd64.tar.gz

        # Configure Go environment variables and apply them
        echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
        echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
        source ~/.bashrc
        ```

      * **Step 3: Download, Build, and Install Singularity**

        ```bash
        # Navigate to a suitable directory for downloading source code
        cd /tmp

        # Download the Singularity CE source code (check for the latest version)
        wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz

        # Extract the archive and clean up
        tar -xvzf singularity-ce-4.0.1.tar.gz
        rm singularity-ce-4.0.1.tar.gz
        cd singularity-ce-4.0.1

        # Configure, build, and install Singularity
        ./mconfig
        cd builddir
        make
        sudo make install
        ```

      * **Step 4: Verify the Installation**

        ```bash
        # Check the installed version
        singularity --version
        ```

3.  **Pipeline Files**:

      * `run_pipeline.sh`
      * `csRNA.sif` (The Singularity container)

4.  **Reference Data**: A directory containing all necessary reference files.

## Setup

### 1\. Prepare the Sample Sheet

This is the most critical input file. Create a CSV file named `samplesheet.csv`. The structure is designed to link csRNA, input, and RNA-seq files for each biological replicate and supports both single-end and paired-end RNA-seq (paried RNA-seq is optional).

* `sample`: A unique identifier for the biological replicate (e.g., `Control_Rep1`).
* `condition`: The experimental group (e.g., `Control`, `Treated`). Used for differential expression.
* `type`: The type of data. Must be one of `csRNA`, `input`, or `rna`.
* `fastq_r1`: The **absolute path** to the R1 FASTQ file.
* `fastq_r2`: The **absolute path** to the R2 FASTQ file. **Leave this column empty for single-end data.**
* `strandedness`: **Required only for `rna` type**. Specifies the library strandedness for `makeTagDirectory`. Options: `SE_plus`, `PE_plus`, `SE_minus`, `PE_minus`.

**Example `samplesheet.csv`:**

```csv
sample,condition,type,fastq_r1,fastq_r2,strandedness
Control_Rep1,Control,csRNA,/path/to/C1_csRNA.fastq.gz,,
Control_Rep1,Control,input,/path/to/C1_input.fastq.gz,,
Control_Rep1,Control,rna,/path/to/C1_rna_R1.fastq.gz,/path/to/C1_rna_R2.fastq.gz,PE_plus
Control_Rep2,Control,csRNA,/path/to/C2_csRNA.fastq.gz,,
Control_Rep2,Control,input,/path/to/C2_input.fastq.gz,,
Treated_Rep1,Treated,csRNA,/path/to/T1_csRNA.fastq.gz,,
Treated_Rep1,Treated,input,/path/to/T1_input.fastq.gz,,
Treated_Rep2,Treated,csRNA,/path/to/T1_csRNA.fastq.gz,,
Treated_Rep2,Treated,input,/path/to/T1_input.fastq.gz,,
Treated_Rep1,Treated,rna,/path/to/T1_rna.fastq.gz,/path/to/C1_rna_T1.fastq.gz,PE_plus
````

### 2\. Prepare the Reference Data

The pipeline requires a genome FASTA, an rRNA FASTA, a gene annotation GTF file, and pre-built STAR indices for both the genome and rRNA sequences.

#### Create Reference Directory & Download Files

```bash
mkdir -p reference_data
cd reference_data

# 1. Download Genome FASTA (from GENCODE)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gzprimary_assembly.genome.fa.gz)
gunzip GRCh38.primary_assembly.genome.fa.gz

# 2. Download Gene Annotation GTF (from GENCODE)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz
gunzip gencode.v46.primary_assembly.annotation.gtf.gz

# 3. Download rRNA sequences (Example from Ensembl)
# (Note: Optional, only for RNA-seq provided)
wget http://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
# Create a list of rRNA intervals from the GTF for bowtie2's rRNA index
awk '$3=="gene" && (/gene_type "rRNA_pseudogene"/ || /gene_type "rRNA"/) {print $1":"$4"-"$5}' gencode.v46.primary_assembly.annotation.gtf > gencode.v46.primary_assembly.rRNA.annotation.txt
singularity exec csRNA.sif samtools faidx GRCh38.primary_assembly.genome.fa -r gencode.v48.primary_assembly.rRNA.annotation.txt > GRCh38.primary_assembly.rRNA.fa
```

#### Build STAR Indices

These are one-time, resource-intensive steps.

```bash
# 1. Build Genome Index
singularity exec ../csRNA.sif STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir ./STAR_index \
  --genomeFastaFiles ./GRCh38.primary_assembly.genome.fa \
  --sjdbGTFfile ./gencode.v46.primary_assembly.annotation.gtf \
  --sjdbOverhang 149

# 2. Build rRNA Index
# (Note: Optional, only for RNA-seq provided)
singularity exec ../csRNA.sif STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir ./rRNA_STAR_index \
  --genomeFastaFiles GRCh38.primary_assembly.rRNA.fa
```

## Running

Execute the pipeline with a single command.

```bash
bash csRNA.sh \
  -s ./samplesheet.csv \
  -o ./csrna_homer_results \
  -r ./reference_data \
  -c ./csRNA.sif \
  -t 8
```

## Output Structure and Interpretation

The output directory will contain results organized by analysis type.

```
./csrna_homer_results/
├── tag_directories/
│   ├── Control_Rep1_csRNA/
│   ├── Control_Rep1_input/
│   └── ...
├── tss_calling/
│   ├── Control_Rep1.tss.txt
│   └── ...
├── differential_expression/
│   ├── merged_tss_clusters.txt
│   ├── raw_counts.txt
│   └── differential_results.txt
├── motif_analysis/
│   ├── nucleotide_freq.txt
│   ├── knownResults.html
│   ├── homerResults.html
│   └── ...
├── visualization/
│   ├── Control_Rep1_csRNA.pos.bedGraph
│   └── Control_Rep1_csRNA.neg.bedGraph
    └── ...
└── multiqc_report/
    └── multiqc_report.html
```
### 1. Quality Control and Alignment

* **`multiqc_report/multiqc_report.html`**: Open this file in a web browser to explore all quality control sections interactively.

    * **Application**: This is the first file you should check to assess the overall quality of your sequencing data and the alignment process. It helps identify problematic samples (e.g., low alignment rate, high duplication) early on.

        * **General Statistics**: A combined table summarizing important metrics for each sample, including read counts and alignment percentages.

        <img width="1986" height="1196" alt="CleanShot 2025-09-25 at 18 35 46@2x" src="https://github.com/user-attachments/assets/f56c93cb-48b0-47d2-ae21-a2f86d842572" />


        * **FastQC**: Reports quality-control metrics on raw and trimmed reads, including 'Sequence Counts', 'Sequence Quality Histograms', 'Per Sequence Quality Scores', and 'Adapter Content'.
        *   - **Sequence Quality Histograms**: The mean quality value across each base position in the read.
         
              <img width="2000" height="1252" alt="CleanShot 2025-09-25 at 18 33 24@2x" src="https://github.com/user-attachments/assets/21bdc1ca-6c77-4681-ba6a-fdd3247a56b9" />


            - **Adapter Content**: The cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position.
         
              <img width="1996" height="1236" alt="CleanShot 2025-09-25 at 18 34 03@2x" src="https://github.com/user-attachments/assets/2ea454a6-8662-4ca1-9b61-7ddcd4b4cb2a" />

        - **Cutadapt**: Reports the number of reads and bases trimmed for adapters and quality:

        <img width="2010" height="1232" alt="CleanShot 2025-09-25 at 18 37 02@2x" src="https://github.com/user-attachments/assets/9c19b496-9ece-47db-88c7-bacb0bd36f46" />

     
          
        * **STAR**: Summarizes alignment metrics, including the percentage of uniquely mapped reads, multimapped reads, and unmapped reads. For the RNA-seq samples, this also reports the percentage of reads removed during rRNA filtering.
     
        * <img width="2006" height="838" alt="CleanShot 2025-09-25 at 18 35 23@2x" src="https://github.com/user-attachments/assets/d9f8ed15-5f3c-4a01-b6d9-b9f6f2a7ba16" />

     
        

* **`tag_directories/`**: These directories are the primary data format used by **HOMER** for all downstream analysis.

    * **Application**: These directories store aligned read positions in an indexed, compressed format, along with crucial metadata (e.g., total mapped reads, fragment length) used by **HOMER** for normalization and peak/TSS calling.

***

### 2. TSS Calling and Quantification

* **`tss_calling/*.tss.txt`**: A Tab-Separated Value file containing the TSS clusters identified for a single sample (e.g., `Control_Rep1.tss.txt`).

    * **Format**: A standard HOMER peak file format. Key columns include:
        * `PeakID`: The ID assigned by `findcsRNATSS.pl`.
        * `Chr`, `Start`, `End`, `Strand`: The genomic location of the TSS cluster.
        * `Annotation`: Genomic feature the TSS is near (e.g., Gene, TSS).
        * `Score`: A statistical measure of significance for the TSS call.
        * `Total Tags`: The raw read count within the cluster for that specific sample.

    * **Application**: This file represents the raw, per-sample output of the **TSS calling** step (`findcsRNATSS.pl`). It is useful for inspecting TSSs specific to a single biological replicate before merging.

* **`differential_expression/merged_tss_clusters.txt`**: A single HOMER peak file containing the combined coordinates of all TSS clusters identified across all samples.

    * **Application**: This file serves as the **master coordinate set** for differential analysis. All subsequent quantification steps use these regions to count reads from all samples.

* **`differential_expression/raw_counts.txt`**: A Tab-Separated Value file generated by `annotatePeaks.pl`, which contains the annotation of the merged TSS clusters and the raw read counts for *every* csRNA sample.

    * **Application**: This is the direct input matrix for the **DESeq2** R analysis. It ensures that all samples have read counts quantified across the identical set of merged genomic regions.

***

### 3. Differential Expression Analysis

* **`differential_expression/differential_results.txt`**: The final output file containing the results of the **DESeq2** analysis.

    * **Format**: A wide, Tab-Separated Value file. It combines the original annotation and raw count columns from `raw_counts.txt` with the statistical results from **DESeq2**. Key columns (from DESeq2) include:
        * `baseMean`: The mean of normalized counts across all samples.
        * `log2FoldChange`: The estimated log2 fold change between the two conditions (e.g., Treated vs. Control).
        * `lfcSE`: The standard error of the log2 fold change.
        * `stat`: The Wald test statistic.
        * `pvalue`: The nominal p-value.
        * **`padj`**: The adjusted p-value (FDR or $q$-value), corrected for multiple testing.

    * **Application**: This is the **primary result file** for the entire pipeline. You will filter this file by `padj` (typically $<0.05$) and `log2FoldChange` (e.g., $|\text{log2FoldChange}| > 1$) to identify **differentially regulated TSS clusters** between your experimental conditions.

***

### 4. Visualization and Motif Analysis

* **`visualization/*.bedGraph`**: Strand-specific **BedGraph** files (e.g., `Control_Rep1_csRNA.pos.bedGraph`).

    * **Format**: A four-column, space-separated file suitable for visualization: `chromosome`, `start`, `end`, `value` (read count). One file is generated for the positive strand (`.pos.bedGraph`) and one for the negative strand (`.neg.bedGraph`) for each sample.

    * **Application**: These files allow you to load the strand-specific csRNA signal for individual samples directly into a **genome browser** (like IGV or UCSC Genome Browser). This is crucial for visually confirming the TSS clusters and observing signal changes between conditions.

    <img width="1968" height="408" alt="CleanShot 2025-09-25 at 18 26 32@2x" src="https://github.com/user-attachments/assets/1b860be0-28ae-4f06-b286-49469c840cc2" />


* **`motif_analysis/knownResults.html`**: 

    * **Application**: This report is used for **motif enrichment analysis**. It provides a ranked list of known transcription factor (TF) binding sites and *de novo* motifs that are statistically enriched near the **merged TSS clusters** compared to a background set of sequences. This helps predict the TFs that may be driving the differential gene expression.
 
    * <img width="2654" height="618" alt="CleanShot 2025-09-25 at 18 27 38@2x" src="https://github.com/user-attachments/assets/dcab98e5-a5b4-4154-affc-40509a6410ff" />
    

* **`motif_analysis/homerResults.html`**: 

    * **Application**: This report is used for **motif enrichment analysis**. It provides a ranked list of homer detected transcription factor (TF) binding sites and *de novo* motifs that are statistically enriched near the **merged TSS clusters** compared to a background set of sequences. This helps predict the TFs that may be driving the differential gene expression.

    <img width="2644" height="784" alt="CleanShot 2025-09-25 at 18 29 15@2x" src="https://github.com/user-attachments/assets/2dd6ccad-6996-4579-902b-8e214564aef4" />
    

* **`motif_analysis/nucleotide_freq.txt`**: A text file containing the nucleotide (A, C, G, T) frequency distribution around the merged TSS clusters.

    * **Application**: Useful for characterizing the general sequence properties of the csRNA-seq peaks, such as the presence of a C/G-rich region or specific initiator elements around the TSS. 
