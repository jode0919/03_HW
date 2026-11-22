import re
import pandas as pd
from snakemake.exceptions import WorkflowError

OUT = config["outdir"].rstrip("/")
REF = config["reference"]
THR = config["threads"]
RES = config["resources"]

# ---------------------------------------------------------------- sample sheet
samples = pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#").fillna("")
samples = samples.apply(lambda col: col.str.strip())
for col in ("sample", "fastq1"):
    if col not in samples.columns:
        raise WorkflowError(f"Sample sheet {config['samples']} needs a '{col}' column")
if "fastq2" not in samples.columns:
    samples["fastq2"] = ""
if samples["sample"].duplicated().any():
    raise WorkflowError("Duplicate sample names in the sample sheet")
if (samples["sample"] == "").any() or (samples["fastq1"] == "").any():
    raise WorkflowError("Every sample needs a name and a fastq1 file")
samples = samples.set_index("sample", drop=False)

SAMPLES = list(samples.index)
PAIRED = [s for s in SAMPLES if samples.loc[s, "fastq2"] != ""]
SINGLE = [s for s in SAMPLES if s not in PAIRED]


def sample_constraint(names):
    return "|".join(re.escape(s) for s in names) or "__none__"


# ------------------------------------------------------- shared path templates
BAM = OUT + "/star/{sample}/Aligned.sortedByCoord.out.bam"
BAI = BAM + ".bai"
BED12 = REF["bed12"]
GENEBODY_BED = config["rseqc"]["genebody_bed"] or BED12
STAR_SA = REF["star_index"] + "/SA"


def mates(sample):
    return ["R1", "R2"] if sample in PAIRED else ["R1"]


RAW_PAIRS = [(s, m) for s in SAMPLES for m in mates(s)]
RAW_S = [p[0] for p in RAW_PAIRS]
RAW_M = [p[1] for p in RAW_PAIRS]


def raw_fastq(wildcards):
    col = "fastq1" if wildcards.mate == "R1" else "fastq2"
    return samples.loc[wildcards.sample, col]


def trimmed_reads(wildcards):
    base = f"{OUT}/trimmed/{wildcards.sample}"
    reads = [base + "_R1.fastq.gz"]
    if wildcards.sample in PAIRED:
        reads.append(base + "_R2.fastq.gz")
    return reads


# ----------------------------------------------------------- target collections
def multiqc_inputs(wildcards=None):
    f = []
    f += expand(OUT + "/fastqc/raw/{sample}_{mate}_fastqc.zip", zip, sample=RAW_S, mate=RAW_M)
    f += expand(OUT + "/fastqc/trimmed/{sample}_{mate}_fastqc.zip", zip, sample=RAW_S, mate=RAW_M)
    f += expand(OUT + "/fastp/{sample}.json", sample=SAMPLES)
    f += expand(OUT + "/star/{sample}/Log.final.out", sample=SAMPLES)
    f += expand(OUT + "/featurecounts/{sample}.counts.txt.summary", sample=SAMPLES)
    f += expand(OUT + "/rseqc/bam_stat/{sample}.txt", sample=SAMPLES)
    f += expand(OUT + "/rseqc/infer_experiment/{sample}.txt", sample=SAMPLES)
    f += expand(OUT + "/rseqc/read_distribution/{sample}.txt", sample=SAMPLES)
    f += expand(OUT + "/rseqc/genebody_coverage/{sample}.geneBodyCoverage.txt", sample=SAMPLES)
    f += expand(OUT + "/rseqc/inner_distance/{sample}.inner_distance_freq.txt", sample=PAIRED)
    f += expand(OUT + "/deeptools/fragment_size/{sample}.table.tsv", sample=PAIRED)
    f += expand(OUT + "/kallisto/{sample}/run_info.json", sample=SAMPLES)
    f += [OUT + "/deeptools/correlation_matrix.tsv", OUT + "/deeptools/pca_data.tsv"]
    return f


def final_targets():
    t = [
        OUT + "/multiqc/multiqc_report.html",
        OUT + "/featurecounts/counts_matrix.tsv",
        OUT + "/kallisto/tpm_matrix.tsv",
        OUT + "/kallisto/est_counts_matrix.tsv",
        OUT + "/deeptools/correlation_heatmap.png",
        OUT + "/deeptools/pca.png",
    ]
    t += expand(OUT + "/rseqc/genebody_coverage/{sample}.geneBodyCoverage.curves.pdf", sample=SAMPLES)
    t += expand(OUT + "/rseqc/inner_distance/{sample}.inner_distance_plot.pdf", sample=PAIRED)
    t += expand(OUT + "/deeptools/fragment_size/{sample}.png", sample=PAIRED)
    if config["deeptools"]["bigwig"]:
        t += expand(OUT + "/deeptools/bigwig/{sample}.bw", sample=SAMPLES)
    return t
