# 2.3 Alignment-free quantification with kallisto


def kallisto_mode(wildcards):
    if wildcards.sample in PAIRED:
        return ""
    k = config["kallisto"]
    return f"--single -l {k['fragment_length_mean']} -s {k['fragment_length_sd']}"


rule kallisto_index:
    """Built only if the index is missing; shared by all samples."""
    input:
        REF["transcripts_fasta"],
    output:
        REF["kallisto_index"],
    log:
        OUT + "/logs/kallisto_index/kallisto_index.log",
    threads: THR["kallisto_index"]
    resources:
        mem_mb=16000,
    conda:
        "../envs/kallisto.yaml"
    shell:
        "kallisto index -i {output} {input} > {log} 2>&1"


rule kallisto_quant:
    input:
        index=REF["kallisto_index"],
        reads=trimmed_reads,
        strand=OUT + "/rseqc/strandedness/{sample}.txt",
    output:
        tsv=OUT + "/kallisto/{sample}/abundance.tsv",
        h5=OUT + "/kallisto/{sample}/abundance.h5",
        info=OUT + "/kallisto/{sample}/run_info.json",
    params:
        outdir=OUT + "/kallisto/{sample}",
        boot=config["kallisto"]["bootstraps"],
        mode=kallisto_mode,
        extra=config["kallisto"]["extra"],
    log:
        OUT + "/logs/kallisto/{sample}.log",
    threads: THR["kallisto"]
    resources:
        mem_mb=RES["default_mem_mb"] * 2,
    conda:
        "../envs/kallisto.yaml"
    shell:
        """
        case "$(cat {input.strand})" in
            forward) ST="--fr-stranded";;
            reverse) ST="--rf-stranded";;
            *) ST="";;
        esac
        kallisto quant -i {input.index} -o {params.outdir} -t {threads} \
                       -b {params.boot} $ST {params.mode} {params.extra} \
                       {input.reads} 2> {log}
        """


rule merge_kallisto:
    input:
        expand(OUT + "/kallisto/{sample}/abundance.tsv", sample=SAMPLES),
    output:
        tpm=OUT + "/kallisto/tpm_matrix.tsv",
        counts=OUT + "/kallisto/est_counts_matrix.tsv",
    params:
        samples=SAMPLES,
    log:
        OUT + "/logs/merge_kallisto/merge_kallisto.log",
    resources:
        mem_mb=RES["default_mem_mb"],
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/merge_kallisto.py \
            --samples {params.samples} --inputs {input} \
            --tpm {output.tpm} --counts {output.counts} > {log} 2>&1
        """
