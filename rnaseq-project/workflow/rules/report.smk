# Aggregate report over all steps


rule multiqc:
    input:
        multiqc_inputs,
    output:
        OUT + "/multiqc/multiqc_report.html",
    params:
        outdir=OUT + "/multiqc",
        dirs=" ".join(
            [
                OUT + "/fastqc",
                OUT + "/fastp",
                OUT + "/star",
                OUT + "/featurecounts",
                OUT + "/rseqc",
                OUT + "/deeptools",
                OUT + "/logs/kallisto",
            ]
        ),
    log:
        OUT + "/logs/multiqc/multiqc.log",
    threads: THR["multiqc"]
    resources:
        mem_mb=RES["default_mem_mb"],
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc --force --outdir {params.outdir} --filename multiqc_report {params.dirs} > {log} 2>&1"
