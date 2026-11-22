RNA-seq Snakemake workflow (human, GRCh38 / GENCODE)
====================================================

Datasets (GEO): GSM1275870, GSM1275866, GSM1275871, GSM1275872

WHAT IT DOES
  2.1 Preprocessing   FastQC (raw) -> fastp (adapter + quality trimming) -> FastQC (trimmed)
  2.2 Alignment       STAR index (built once if missing) -> STAR -> featureCounts (+ count matrix)
  2.3 Alignment-free  kallisto index (built once if missing) -> kallisto quant (+ TPM/count matrices)
  2.4 Alignment QC    RSeQC: bam_stat, infer_experiment (strandedness), read_distribution,
                      geneBody_coverage, inner_distance
                      deepTools: bamPEFragmentSize, multiBamSummary -> plotCorrelation, plotPCA, bigWig
      Everything is summarised in results/multiqc/multiqc_report.html

LAYOUT
  workflow/Snakefile, workflow/rules/*.smk   Snakemake rules
  workflow/envs/*.yaml                       one Conda environment per tool
  workflow/scripts/                          small Python helpers used by rules
  config/config.yaml                         parameters, reference paths, threads, memory
  config/samples.tsv                         sample sheet (sample, fastq1, fastq2)
  scripts/download_sra.sh                    GSM -> SRR -> FASTQ (SRA Toolkit)
  scripts/prepare_reference.sh               GENCODE genome/annotation/transcriptome download
  scripts/make_samples_tsv.sh                build samples.tsv from a FASTQ directory

REQUIREMENTS
  Conda/Miniforge with Snakemake (>= 8):  conda create -n snakemake -c conda-forge -c bioconda snakemake
  SRA Toolkit for download_sra.sh (conda install -c bioconda sra-tools)
  STAR human index needs ~32-40 GB RAM; alignment ~30 GB. Disk: >= 100 GB recommended.

HOW TO RUN
  1. bash scripts/download_sra.sh accessions.txt        # accessions.txt = GSM ids, one per line
  2. bash scripts/prepare_reference.sh                  # or: GENCODE_RELEASE=47 bash ...
  3. bash scripts/make_samples_tsv.sh fastq > config/samples.tsv
  4. snakemake --use-conda --cores 16 --resources mem_mb=48000 -n     # dry run
     snakemake --use-conda --cores 16 --resources mem_mb=48000        # real run
  New dataset (same species/style): change config/samples.tsv (and config.yaml if needed);
  the Snakefile is not modified. Single-end and paired-end samples are detected automatically
  (empty fastq2 = single-end); paired-only analyses run only for paired samples.
  Strandedness is inferred per sample (config: strandedness: auto) or can be forced.
  Logs: results/logs/<rule>/<sample>.log

WHY THESE GENCODE FILES (GRCh38)
  Genome  GRCh38.primary_assembly.genome.fa
          Primary assembly = chromosomes + unplaced scaffolds, WITHOUT ALT/patch haplotypes.
          ALT contigs duplicate sequence and cause artificial multi-mapping; the "ALL"
          file is therefore avoided. Unlike the chromosome-only file it keeps unplaced scaffolds
          so reads from those loci are not forced onto wrong chromosomes.
  GTF     gencode.vNN.primary_assembly.annotation.gtf  (comprehensive, PRI regions)
          Same coordinates and chromosome names as the genome file; the comprehensive set
          includes all transcripts (lncRNA, etc.), unlike the smaller "basic" set. Used for
          STAR splice junctions, featureCounts, and the RSeQC BED12 model.
  Transcriptome  gencode.vNN.transcripts.fa.gz
          Transcript sequences from the same release/annotation, used for kallisto so both
          quantification routes refer to identical gene models.
  Pin the release with GENCODE_RELEASE=NN for reproducibility; VERSION.txt records what was used.

NOT IN GIT
  FASTQ/SRA files, reference files/indices and results are excluded via .gitignore.
