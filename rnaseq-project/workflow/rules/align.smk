# 2.2 Alignment-based quantification: STAR index (once) -> STAR -> featureCounts


rule star_index:
    """Built only if the index is missing; shared by all samples."""
    input:
        fasta=REF["genome_fasta"],
        gtf=REF["gtf"],
    output:
        sa=STAR_SA,
        genome=REF["star_index"] + "/Genome",
        sai=REF["star_index"] + "/SAindex",
    params:
        dir=REF["star_index"],
        overhang=config["star"]["sjdb_overhang"],
        extra=config["star"].get("index_extra", ""),
        ram=lambda wc, resources: int(resources.mem_mb * 0.8) * 1_000_000,
    log:
        OUT + "/logs/star_index/star_index.log",
    threads: THR["star_index"]
    resources:
        mem_mb=RES["star_index_mem_mb"],
    conda:
        "../envs/star.yaml"
    shell:
        """
        STAR --runMode genomeGenerate --runThreadN {threads} \
             --genomeDir {params.dir} --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} --sjdbOverhang {params.overhang} \
             --limitGenomeGenerateRAM {params.ram} {params.extra} \
             --outFileNamePrefix {params.dir}/ > {log} 2>&1
        """


rule star_align:
    input:
        reads=trimmed_reads,
        index=STAR_SA,
    output:
        bam=BAM,
        final=OUT + "/star/{sample}/Log.final.out",
        sj=OUT + "/star/{sample}/SJ.out.tab",
    params:
        index=REF["star_index"],
        prefix=OUT + "/star/{sample}/",
        tmp=OUT + "/star/{sample}_STARtmp",
        extra=config["star"]["extra"],
    log:
        OUT + "/logs/star/{sample}.log",
    threads: THR["star"]
    resources:
        mem_mb=RES["star_mem_mb"],
    conda:
        "../envs/star.yaml"
    shell:
        """
        rm -rf {params.tmp}
        STAR --runThreadN {threads} --genomeDir {params.index} \
             --readFilesIn {input.reads} --readFilesCommand zcat \
             --outFileNamePrefix {params.prefix} --outTmpDir {params.tmp} \
             --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM \
             {params.extra} > {log} 2>&1
        """


rule samtools_index:
    input:
        BAM,
    output:
        BAI,
    log:
        OUT + "/logs/samtools_index/{sample}.log",
    threads: THR["samtools"]
    resources:
        mem_mb=RES["default_mem_mb"],
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -@ {threads} {input} {output} 2> {log}"


rule featurecounts:
    """Per-sample counts; strandedness (-s) comes from RSeQC / config."""
    input:
        bam=BAM,
        bai=BAI,
        gtf=REF["gtf"],
        strand=OUT + "/rseqc/strandedness/{sample}.txt",
    output:
        counts=OUT + "/featurecounts/{sample}.counts.txt",
        summary=OUT + "/featurecounts/{sample}.counts.txt.summary",
    params:
        paired=lambda wc: "-p --countReadPairs" if wc.sample in PAIRED else "",
        ft=config["featurecounts"]["feature_type"],
        attr=config["featurecounts"]["attribute"],
        extra=config["featurecounts"]["extra"],
    log:
        OUT + "/logs/featurecounts/{sample}.log",
    threads: THR["featurecounts"]
    resources:
        mem_mb=RES["default_mem_mb"],
    conda:
        "../envs/subread.yaml"
    shell:
        """
        case "$(cat {input.strand})" in forward) S=1;; reverse) S=2;; *) S=0;; esac
        featureCounts -T {threads} -a {input.gtf} -o {output.counts} \
                      -t {params.ft} -g {params.attr} -s $S {params.paired} \
                      {params.extra} {input.bam} > {log} 2>&1
        """


rule merge_counts:
    input:
        expand(OUT + "/featurecounts/{sample}.counts.txt", sample=SAMPLES),
    output:
        OUT + "/featurecounts/counts_matrix.tsv",
    params:
        samples=SAMPLES,
    log:
        OUT + "/logs/merge_counts/merge_counts.log",
    resources:
        mem_mb=RES["default_mem_mb"],
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/merge_counts.py \
            --samples {params.samples} --inputs {input} --out {output} > {log} 2>&1
        """
