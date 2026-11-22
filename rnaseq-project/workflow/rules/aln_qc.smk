# 2.3 Alignment quality assessment: RSeQC + deepTools


rule gtf_to_bed12:
    """RSeQC needs a BED12 gene model; built from the GTF only if missing."""
    input:
        REF["gtf"],
    output:
        BED12,
    log:
        OUT + "/logs/gtf_to_bed12/gtf_to_bed12.log",
    resources:
        mem_mb=RES["default_mem_mb"],
    conda:
        "../envs/ucsc.yaml"
    shell:
        """
        tmp=$(mktemp)
        gtfToGenePred -genePredExt -geneNameAsName2 {input} $tmp 2> {log}
        genePredToBed $tmp {output} 2>> {log}
        rm -f $tmp
        """


# ------------------------------------------------------------------- RSeQC
rule rseqc_bam_stat:
    input:
        bam=BAM,
        bai=BAI,
    output:
        OUT + "/rseqc/bam_stat/{sample}.txt",
    log:
        OUT + "/logs/rseqc_bam_stat/{sample}.log",
    threads: THR["rseqc"]
    resources:
        mem_mb=RES["default_mem_mb"],
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input.bam} > {output} 2> {log}"


rule rseqc_infer_experiment:
    input:
        bam=BAM,
        bai=BAI,
        bed=BED12,
    output:
        OUT + "/rseqc/infer_experiment/{sample}.txt",
    params:
        n=config["rseqc"]["infer_sample_size"],
    log:
        OUT + "/logs/rseqc_infer_experiment/{sample}.log",
    threads: THR["rseqc"]
    resources:
        mem_mb=RES["default_mem_mb"],
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -i {input.bam} -r {input.bed} -s {params.n} > {output} 2> {log}"


rule call_strandedness:
    """Turn infer_experiment output into unstranded/forward/reverse (or use the config value)."""
    input:
        OUT + "/rseqc/infer_experiment/{sample}.txt",
    output:
        OUT + "/rseqc/strandedness/{sample}.txt",
    params:
        mode=config["strandedness"],
        thr=config["strand_threshold"],
    log:
        OUT + "/logs/call_strandedness/{sample}.log",
    resources:
        mem_mb=500,
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/call_strand.py --infer {input} \
            --mode {params.mode} --threshold {params.thr} --out {output} > {log} 2>&1
        """


rule rseqc_read_distribution:
    input:
        bam=BAM,
        bai=BAI,
        bed=BED12,
    output:
        OUT + "/rseqc/read_distribution/{sample}.txt",
    log:
        OUT + "/logs/rseqc_read_distribution/{sample}.log",
    threads: THR["rseqc"]
    resources:
        mem_mb=RES["default_mem_mb"] * 2,
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -i {input.bam} -r {input.bed} > {output} 2> {log}"


rule rseqc_genebody_coverage:
    input:
        bam=BAM,
        bai=BAI,
        bed=GENEBODY_BED,
    output:
        txt=OUT + "/rseqc/genebody_coverage/{sample}.geneBodyCoverage.txt",
        pdf=OUT + "/rseqc/genebody_coverage/{sample}.geneBodyCoverage.curves.pdf",
    params:
        prefix=OUT + "/rseqc/genebody_coverage/{sample}",
    log:
        OUT + "/logs/rseqc_genebody_coverage/{sample}.log",
    threads: THR["rseqc"]
    resources:
        mem_mb=RES["default_mem_mb"] * 2,
    conda:
        "../envs/rseqc.yaml"
    shell:
        "geneBody_coverage.py -i {input.bam} -r {input.bed} -o {params.prefix} > {log} 2>&1"


rule rseqc_inner_distance:
    """Paired-end only."""
    input:
        bam=BAM,
        bai=BAI,
        bed=BED12,
    output:
        freq=OUT + "/rseqc/inner_distance/{sample}.inner_distance_freq.txt",
        pdf=OUT + "/rseqc/inner_distance/{sample}.inner_distance_plot.pdf",
    wildcard_constraints:
        sample=sample_constraint(PAIRED),
    params:
        prefix=OUT + "/rseqc/inner_distance/{sample}",
        k=config["rseqc"]["inner_distance_reads"],
    log:
        OUT + "/logs/rseqc_inner_distance/{sample}.log",
    threads: THR["rseqc"]
    resources:
        mem_mb=RES["default_mem_mb"],
    conda:
        "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -i {input.bam} -r {input.bed} -k {params.k} -o {params.prefix} > {log} 2>&1"


# ---------------------------------------------------------------- deepTools
rule deeptools_fragment_size:
    """Fragment length distribution, paired-end only."""
    input:
        bam=BAM,
        bai=BAI,
    output:
        png=OUT + "/deeptools/fragment_size/{sample}.png",
        table=OUT + "/deeptools/fragment_size/{sample}.table.tsv",
    wildcard_constraints:
        sample=sample_constraint(PAIRED),
    log:
        OUT + "/logs/deeptools_fragment_size/{sample}.log",
    threads: THR["deeptools"]
    resources:
        mem_mb=RES["deeptools_mem_mb"],
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        bamPEFragmentSize -b {input.bam} --samplesLabel {wildcards.sample} \
            -hist {output.png} -T "Fragment size: {wildcards.sample}" \
            --table {output.table} -p {threads} > {log} 2>&1
        """


rule deeptools_multibamsummary:
    input:
        bams=expand(BAM, sample=SAMPLES),
        bais=expand(BAI, sample=SAMPLES),
    output:
        npz=OUT + "/deeptools/multiBamSummary.npz",
        counts=OUT + "/deeptools/multiBamSummary.counts.tsv",
    params:
        labels=" ".join(SAMPLES),
        bs=config["deeptools"]["bin_size"],
    log:
        OUT + "/logs/deeptools_multibamsummary/multibamsummary.log",
    threads: THR["deeptools"]
    resources:
        mem_mb=RES["deeptools_mem_mb"],
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        multiBamSummary bins -b {input.bams} --labels {params.labels} \
            -bs {params.bs} -p {threads} -o {output.npz} \
            --outRawCounts {output.counts} > {log} 2>&1
        """


rule deeptools_plot_correlation:
    input:
        OUT + "/deeptools/multiBamSummary.npz",
    output:
        png=OUT + "/deeptools/correlation_heatmap.png",
        matrix=OUT + "/deeptools/correlation_matrix.tsv",
    params:
        method=config["deeptools"]["cor_method"],
    log:
        OUT + "/logs/deeptools_plot_correlation/plot_correlation.log",
    resources:
        mem_mb=RES["deeptools_mem_mb"],
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        plotCorrelation -in {input} --corMethod {params.method} --skipZeros \
            --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
            -o {output.png} --outFileCorMatrix {output.matrix} > {log} 2>&1
        """


rule deeptools_plot_pca:
    input:
        OUT + "/deeptools/multiBamSummary.npz",
    output:
        png=OUT + "/deeptools/pca.png",
        data=OUT + "/deeptools/pca_data.tsv",
    log:
        OUT + "/logs/deeptools_plot_pca/plot_pca.log",
    resources:
        mem_mb=RES["deeptools_mem_mb"],
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotPCA -in {input} -o {output.png} --outFileNameData {output.data} > {log} 2>&1"


rule deeptools_bigwig:
    input:
        bam=BAM,
        bai=BAI,
    output:
        OUT + "/deeptools/bigwig/{sample}.bw",
    params:
        bs=config["deeptools"]["bigwig_bin_size"],
    log:
        OUT + "/logs/deeptools_bigwig/{sample}.log",
    threads: THR["deeptools"]
    resources:
        mem_mb=RES["deeptools_mem_mb"],
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        bamCoverage -b {input.bam} -o {output} --normalizeUsing CPM \
            --binSize {params.bs} -p {threads} > {log} 2>&1
        """
