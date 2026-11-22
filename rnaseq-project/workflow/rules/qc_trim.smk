# 2.1 Read preprocessing: raw QC -> trimming -> QC of trimmed reads


rule fastqc_raw:
    input:
        raw_fastq,
    output:
        html=OUT + "/fastqc/raw/{sample}_{mate}_fastqc.html",
        zip=OUT + "/fastqc/raw/{sample}_{mate}_fastqc.zip",
    wildcard_constraints:
        mate="R[12]",
    log:
        OUT + "/logs/fastqc_raw/{sample}_{mate}.log",
    threads: THR["fastqc"]
    resources:
        mem_mb=RES["default_mem_mb"],
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        tmp=$(mktemp -d)
        fastqc -t {threads} -o $tmp {input} > {log} 2>&1
        mv $tmp/*_fastqc.html {output.html}
        mv $tmp/*_fastqc.zip {output.zip}
        rm -rf $tmp
        """


rule fastp_pe:
    input:
        r1=lambda wc: samples.loc[wc.sample, "fastq1"],
        r2=lambda wc: samples.loc[wc.sample, "fastq2"],
    output:
        r1=OUT + "/trimmed/{sample}_R1.fastq.gz",
        r2=OUT + "/trimmed/{sample}_R2.fastq.gz",
        json=OUT + "/fastp/{sample}.json",
        html=OUT + "/fastp/{sample}.html",
    wildcard_constraints:
        sample=sample_constraint(PAIRED),
    params:
        q=config["trimming"]["min_quality"],
        minlen=config["trimming"]["min_length"],
        extra=config["trimming"]["extra"],
    log:
        OUT + "/logs/fastp/{sample}.log",
    threads: THR["fastp"]
    resources:
        mem_mb=RES["default_mem_mb"],
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
              --json {output.json} --html {output.html} --thread {threads} \
              --detect_adapter_for_pe --cut_tail --cut_mean_quality {params.q} \
              --qualified_quality_phred {params.q} --length_required {params.minlen} \
              {params.extra} 2> {log}
        """


rule fastp_se:
    input:
        r1=lambda wc: samples.loc[wc.sample, "fastq1"],
    output:
        r1=OUT + "/trimmed/{sample}_R1.fastq.gz",
        json=OUT + "/fastp/{sample}.json",
        html=OUT + "/fastp/{sample}.html",
    wildcard_constraints:
        sample=sample_constraint(SINGLE),
    params:
        q=config["trimming"]["min_quality"],
        minlen=config["trimming"]["min_length"],
        extra=config["trimming"]["extra"],
    log:
        OUT + "/logs/fastp/{sample}.log",
    threads: THR["fastp"]
    resources:
        mem_mb=RES["default_mem_mb"],
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        fastp -i {input.r1} -o {output.r1} \
              --json {output.json} --html {output.html} --thread {threads} \
              --cut_tail --cut_mean_quality {params.q} \
              --qualified_quality_phred {params.q} --length_required {params.minlen} \
              {params.extra} 2> {log}
        """


rule fastqc_trimmed:
    input:
        OUT + "/trimmed/{sample}_{mate}.fastq.gz",
    output:
        html=OUT + "/fastqc/trimmed/{sample}_{mate}_fastqc.html",
        zip=OUT + "/fastqc/trimmed/{sample}_{mate}_fastqc.zip",
    wildcard_constraints:
        mate="R[12]",
    log:
        OUT + "/logs/fastqc_trimmed/{sample}_{mate}.log",
    threads: THR["fastqc"]
    resources:
        mem_mb=RES["default_mem_mb"],
    conda:
        "../envs/fastqc.yaml"
    shell:
        """
        tmp=$(mktemp -d)
        fastqc -t {threads} -o $tmp {input} > {log} 2>&1
        mv $tmp/*_fastqc.html {output.html}
        mv $tmp/*_fastqc.zip {output.zip}
        rm -rf $tmp
        """
