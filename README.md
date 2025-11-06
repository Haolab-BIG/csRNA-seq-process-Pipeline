Here is the complete `README.md` file, formatted according to the `README-end-seq_reference.md` template, but using the content and file-specific details from your `csRNA-seq` pipeline, Snakemake file, and R script.

-----

# csRNA-seq-Processing-Pipeline

This csRNA-seq pipeline provides a fully containerized Singularity environment that bundles all required tools (FastQC, Trim Galore, STAR, HOMER, DESeq2). With a single command, the entire workflow—from raw FASTQ input, trimming, quality control, genome alignment, and TSS calling, through to differential expression and motif analysis—can be executed reproducibly on any compatible system.

# Part I Workflow

<img width="2242" height="508" alt="pipeline" src="https://github.com/user-attachments/assets/d39f07eb-e0ab-4499-bc4e-9fee73bfb842" />


# Part II Requirements

1.  **Recommended Specs**:

      * 8-core CPU
      * 64 GB RAM

2.  **Singularity**: Must be installed on your system. Below are the detailed steps for installing on an Ubuntu 22.0.4 system. For other operating systems, please refer to the official installation guide: [https://docs.sylabs.io/guides/latest/user-guide/installation.html](https://www.google.com/search?q=https://docs.sylabs.io/guides/latest/user-guide/installation.html)

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

3.  **snakemake**: Snakemake must be installed on your system and requires a Python 3 distribution.

    ```bash
    pip install snakemake
    ```

4.  **Reference Data**: A directory containing the genome FASTA, GTF, and pre-built STAR indices. Below are detailed steps for the human hg38 genome.

    ```bash
    mkdir reference_data
    cd reference_data

    # 1. Download Genome FASTA (from GENCODE)
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
    gunzip GRCh38.primary_assembly.genome.fa.gz

    # 2. Download Gene Annotation GTF (from GENCODE)
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz
    gunzip gencode.v46.primary_assembly.annotation.gtf.gz

    # 3. Create rRNA FASTA (for optional RNA-seq alignment)
    # (Requires csRNA.sif container to be present in the parent directory)
    awk '$3=="gene" && (/gene_type "rRNA_pseudogene"/ || /gene_type "rRNA"/) {print $1":"$4"-"$5}' gencode.v46.primary_assembly.annotation.gtf > gencode.v46.rRNA.intervals.txt

    singularity exec ../Containers/csRNA.sif samtools faidx GRCh38.primary_assembly.genome.fa -r gencode.v46.rRNA.intervals.txt > GRCh38.primary_assembly.rRNA.fa

    # 4. Build STAR Genome Index (This is resource-intensive)
    mkdir star_genome_index
    singularity exec ../Containers/csRNA.sif STAR \
      --runThreadN 8 \
      --runMode genomeGenerate \
      --genomeDir ./star_genome_index \
      --genomeFastaFiles ./GRCh38.primary_assembly.genome.fa \
      --sjdbGTFfile ./gencode.v46.primary_assembly.annotation.gtf \
      --sjdbOverhang 149

    # 5. Build STAR rRNA Index (for optional RNA-seq alignment)
    mkdir star_rrna_index
    singularity exec ../Containers/csRNA.sif STAR \
      --runThreadN 8 \
      --runMode genomeGenerate \
      --genomeDir ./star_rrna_index \
      --genomeFastaFiles ./GRCh38.primary_assembly.rRNA.fa
    ```

5.  **Required Files**:

    ```bash
    project_directory/
    ├── Scripts/
    │   ├── csRNA-seq.smk
    │   ├── config.yaml
    │   └── csrna_deseq2_analysis.R
    ├── Containers/
    │   └── csRNA.sif
    ├── Reference_Data/
    │   ├── GRCh38.primary_assembly.genome.fa
    │   ├── gencode.v46.primary_assembly.annotation.gtf
    │   ├── GRCh38.primary_assembly.rRNA.fa
    │   ├── star_genome_index/
    │   │   └── (STAR index files)
    │   └── star_rrna_index/
    │       └── (STAR index files)
    └── Raw_Data/
        ├── SRR22752762.fastq.gz
        ├── SRR22752761.fastq.gz
        └── ... (all other fastq files)
    ```

      - **csRNA-seq.smk** — The main Snakemake workflow script.
      - **config.yaml** — Configuration file containing paths, parameters, and sample information.
      - **csrna\_deseq2\_analysis.R** — R script for performing DESeq2 analysis.
      - **csRNA.sif** — Singularity container image with all required software and dependencies pre-installed.
      - **Reference\_Data/** — Directory containing all reference files and STAR indices as built in Step 4.
      - **Raw\_Data/** — Directory containing your raw FASTQ files.

# Part III Running

  * **Example code**

      * **Step 1: Edit `config.yaml`**

        **Note:** All paths must be **absolute paths** and accessible from within the Singularity container (i.e., they must be under the `/project_directory` that you bind-mount).

        ```yaml
        # ============================================================================
        # csRNA-seq Pipeline Configuration File (HOMER-based)
        # ============================================================================

        # --- 1. General Settings ---
        output_dir: "/project_directory/results"

        # Path to the Singularity container (csRNA.sif)
        container: "/project_directory/Containers/csRNA.sif"

        # Number of CPU threads to use for parallel processing
        threads: 8

        # --- 2. Species Configuration ---
        species: "human" # Supported: "human", "mouse"

        species_config:
          human:
            r_annotation_db: "org.Hs.eg.db"

        # --- 3. Reference Files ---
        ref:
          genome_fasta: "/project_directory/Reference_Data/GRCh38.primary_assembly.genome.fa"
          gtf: "/project_directory/Reference_Data/gencode.v46.primary_assembly.annotation.gtf"
          rrna_star_index: "/project_directory/Reference_Data/star_rrna_index"
          star_genome_index: "/project_directory/Reference_Data/star_genome_index"

        # --- 4. Samples ---
        # Define each biological sample.
        # Omit any file type (e.g., 'input' or 'rna') if not available.
        samples:
          Control_Rep1:
            condition: "Control"
            csRNA: "/project_directory/Raw_Data/SRR22752762.fastq.gz"
            input: "/project_directory/Raw_Data/SRR22752761.fastq.gz"
          
          Control_Rep2:
            condition: "Control"
            csRNA: "/project_directory/Raw_Data/SRR22752513.fastq.gz"
            input: "/project_directory/Raw_Data/SRR22752760.fastq.gz"

          Treated_Rep1:
            condition: "Treated"
            csRNA: "/project_directory/Raw_Data/SRR22752511.fastq.gz"
            input: "/project_directory/Raw_Data/SRR22752757.fastq.gz"

          Treated_Rep2:
            condition: "Treated"
            csRNA: "/project_directory/Raw_Data/SRR22752510.fastq.gz"
            input: "/project_directory/Raw_Data/SRR22752756.fastq.gz"

        # --- 5. Analysis Parameters ---
        params:
          trim_galore:
            stringency: 4
            length: 20
            extra_adapter_params: ""

          star_csrna:
            outFilterMultimapNmax: 10000
            outSAMmultNmax: 1

          makeTagDirectory_csrna:
            fragLength: 150

          findMotifsGenome:
            size: "-150,50"
            len: "8,10,12"

        # --- 6. Differential Expression Settings ---
        analysis:
          contrast: ["Treated", "Control"]
          fdr_threshold: 0.05
        ```

      * **Step 2: run snakemake**

        Navigate to the `Scripts/` directory (or wherever `csRNA-seq.smk` and `config.yaml` are located).

        ```bash
        snakemake -s csRNA-seq.smk --use-singularity --cores 8 --singularity-args "--bind /project_directory:/project_directory"
        ```
        Then delete the intermediate files and folders:

        ```bash
        rm -r results/qc/ results/intermediate/
        ```


  * **Command Parameters**

    **edit `config.yaml`**

      - `output_dir`: Path to the directory where all output will be stored (required).
      - `container`: Absolute path to the `csRNA.sif` Singularity container image (required).
      - `threads`: Number of threads to use for multithreaded rules (e.g., STAR, Trim Galore) (required).
      - `species`: Species to use for analysis (e.g., "human", "mouse") (required).
      - `ref`: Section containing absolute paths to all reference files (FASTA, GTF, STAR indices) (required).
      - `samples`: The main sample definition block. Each sample must have a unique name (e.g., `Control_Rep1`) and a `condition` (e.g., `Control`). Provide absolute paths to `csRNA` and `input` fastq files. (Optional `rna` section not shown in example but supported by the workflow).
      - `params`: Section for configuring parameters for specific tools like `trim_galore`, `star_csrna`, and `findMotifsGenome`.
      - `analysis`: Section for defining the differential expression contrast. `contrast: [treatment_group, control_group]` (required for DE analysis).

    **run snakemake**

      - `--use-singularity`: Enables execution of rules within a Singularity container to ensure a fully reproducible environment.
      - `--singularity-args`: Allows passing additional arguments to the Singularity runtime.
      - `--cores`: Specifies the maximum number of CPU cores (threads) that Snakemake can use in parallel.
      - `--bind`: Specifies the directories to be mounted within the Singularity container. **You must bind-mount the entire project directory** so the container can access scripts, raw data, references, and the output directory. The format is `/host/path:/container/path`.

# Part IV. Output

## **Output Structure**

The output directory (defined in `config.yaml`) contains results organized by analysis type:

<output_dir>/
├── tss_calling/
│   ├── Control_Rep1.tss.txt
│   └── … (one for each csRNA sample)
├── results/
│   ├── differential_expression/
│   │   ├── merged_tss_clusters.txt
│   │   ├── raw_counts.txt
│   │   └── deseq2_results.txt
│   ├── motif_analysis/
│   │   ├── nucleotide_freq.txt
│   │   ├── knownResults.html
│   │   ├── homerResults.html
│   │   └── … (motif discovery results)
│   ├── visualization/
│   │   ├── Control_Rep1_csRNA.pos.bedGraph
│   │   ├── Control_Rep1_csRNA.neg.bedGraph
│   │   └── … (strand-specific BedGraphs)
│   └── multiqc_report.html


---

## **Output Interpretation**

### **1. Quality Control and Alignment**

**File:** `multiqc_report/multiqc_report.html`  
**Purpose:** Aggregated quality control summary for raw reads, adapter trimming, and alignment results.

- **General Statistics:** Summary table of sequencing metrics such as total reads, GC content, and alignment rate.  
  ![QC Summary](https://github.com/user-attachments/assets/f56c93cb-48b0-47d2-ae21-a2f86d842572)

- **FastQC:** Per-sample quality profiles.  
  - **Sequence Quality Histograms:** Base-level Phred quality distribution.  
    ![FastQC Quality](https://github.com/user-attachments/assets/21bdc1ca-6c77-4681-ba6a-fdd3247a56b9)
  - **Adapter Content:** Percentage of reads containing adapter sequences.  
    ![Adapter Content](https://github.com/user-attachments/assets/2ea454a6-8662-4ca1-9b61-7ddcd4b4cb2a)

- **Cutadapt:** Adapter and low-quality base trimming statistics.  
  ![Cutadapt](https://github.com/user-attachments/assets/9c19b496-9ece-47db-88c7-bacb0bd36f46)

- **STAR Alignment:** Mapping quality metrics, including unique/multimapped rates and rRNA filtering.  
  ![STAR Alignment](https://github.com/user-attachments/assets/d9f8ed15-5f3c-4a01-b6d9-b9f6f2a7ba16)

---

### **2. TSS Calling and Quantification**

#### **`tss_calling/*.tss.txt`**

- **Description:** Per-sample TSS cluster calls from `findcsRNATSS.pl`.  
- **Key Columns:**
  - `PeakID`: Unique TSS ID  
  - `Chr`, `Start`, `End`, `Strand`: Genomic coordinates  
  - `Annotation`: Gene context (e.g., promoter, exon)  
  - `Score`: Cluster significance metric  
  - `Total Tags`: Raw tag counts per region  

> **Use:** Inspect individual biological replicates before merging for differential testing.

---

#### **`differential_expression/merged_tss_clusters.txt`**

- **Description:** Unified coordinate reference of all TSS clusters across samples.  
- **Use:** Serves as the master feature set for read quantification and DE analysis.

---

#### **`differential_expression/raw_counts.txt`**

- **Description:** Annotated count matrix produced by HOMER’s `annotatePeaks.pl`.  
- **Use:** Direct input for DESeq2 analysis ensuring counts are measured across identical genomic loci in all samples.

---

### **3. Differential Expression Analysis**

#### **`differential_expression/deseq2_results.txt`**

- **Description:** Final differential expression results from DESeq2.  
- **Key Columns:**
  - `baseMean`: Mean normalized counts  
  - `log2FoldChange`: Fold-change between experimental groups  
  - `lfcSE`: Standard error of fold-change  
  - `stat`, `pvalue`, `padj`: Wald test statistic and multiple-testing corrected FDR  

> **Interpretation:**  
> Filter by `padj < 0.05` and `|log2FoldChange| > 1` to identify significantly altered TSS clusters.

---

### **4. Visualization and Motif Analysis**

#### **`visualization/*.bedGraph`**

- **Format:** `chromosome`, `start`, `end`, `value`  
- **Use:** Load `.pos.bedGraph` and `.neg.bedGraph` in IGV or UCSC Genome Browser to visualize strand-specific transcription activity.  

  ![BedGraph](https://github.com/user-attachments/assets/1b860be0-28ae-4f06-b286-49469c840cc2)

---

#### **Motif Enrichment Results**

| File | Description | Use |
|------|--------------|-----|
| `motif_analysis/knownResults.html` | Enrichment of **known TF motifs** near merged TSS clusters. | Identify regulatory TFs contributing to expression changes. |
| `motif_analysis/homerResults.html` | Results of *de novo* motif discovery. | Detect novel or uncharacterized sequence motifs. |
| `motif_analysis/nucleotide_freq.txt` | Base composition around TSS regions. | Examine A/T- or G/C-rich sequence features. |

Example visualizations of **motif enrichment results**:  

![Known Motifs](https://github.com/user-attachments/assets/dcab98e5-a5b4-4154-affc-40509a6410ff)  
![De novo Motifs](https://github.com/user-attachments/assets/2dd6ccad-6996-4579-902b-8e214564aef4)
