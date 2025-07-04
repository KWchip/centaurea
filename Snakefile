configfile: "config/config.yaml"
SAMPLES = config["samples"]
REF     = config["reference"]

rule all:
    input:
        expand("results/variants/{sample}.vcf.gz", sample=SAMPLES)

rule fastqc:
    input:
        r1 = "data/{sample}.1.fq.gz",
        r2 = "data/{sample}.2.fq.gz"
    output:
        html1 = "results/fastqc/{sample}_1.html",
        zip1  = "results/fastqc/{sample}_1_fastqc.zip",
        html2 = "results/fastqc/{sample}_2.html",
        zip2  = "results/fastqc/{sample}_2_fastqc.zip"
    log:
        "logs/fastqc/{sample}.log"
    threads: 1
    conda:
        "envs/fastqc.yaml"
    wrapper:
        "v5.7.0/bio/fastqc"

rule trim_galore_pe:
    input:
        ["data/{sample}.1.fq.gz", "data/{sample}.2.fq.gz"]
    output:
        fasta_fwd = "results/trimmed/{sample}_val_1.fq.gz",
        report_fwd = "results/trimmed/reports/{sample}_R1_trimming_report.txt",
        fasta_rev = "results/trimmed/{sample}_val_2.fq.gz",
        report_rev = "results/trimmed/{sample}_R2_trimming_report.txt"
        # r1_unpaired  = "results/trimmed/{sample}_unpaired_1.fq.gz",
        # r2_unpaired  = "results/trimmed/{sample}_unpaired_2.fq.gz"
    log:
        "logs/trim_galore/{sample}.log"
    threads: 4
    params:
        extra = "--illumina -q 20"
    conda:
        "envs/trim_galore.yaml"
    wrapper:
        "v7.1.0/bio/trim_galore/pe"

rule bwa_index:
    input:
        ref = REF
    output:
        REF + ".amb",
        REF + ".ann",
        REF + ".bwt",
        REF + ".pac",
        REF + ".sa"
    log:
        "logs/bwa/index.log"
    threads: 2
    conda:
        "envs/bwa.yaml"
    shell:
        """
        bwa index -a bwtsw \
            -p {input.ref} \
            {input.ref} \
        2> {log}
        """
        
rule bwa_mem:
    input:
        reads = ["results/trimmed/{sample}_val_1.fq.gz",
                 "results/trimmed/{sample}_val_2.fq.gz"],
        idx    = multiext(REF, ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        "results/aligned/{sample}.bam"
    log:
        "logs/bwa/mem/{sample}.log"
    threads: 8
    params:
        extra      = r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sorting    = "samtools",
        sort_order = "coordinate"
    conda:
        "envs/bwa.yaml"
    wrapper:
        "v1.0.0/bio/bwa/mem" 

rule samtools_faidx:
    input:
        ref = REF
    output:
        REF + ".fai"
    log:
        "logs/samtools/faidx.log"
    threads:1
    conda: 
        "envs/samtools.yaml"
    shell:
        """
        samtools faidx {input.ref} \
        2> {log}
        """

rule bcftools_mpileup:
    input:
        index = REF + ".fai",
        ref   = REF,
        alignments   = "results/aligned/{sample}.bam"
    output:
        pileup = "results/variants/{sample}.bcf"
    log:
        "logs/bcftools/mpileup/{sample}.log"
    threads: 4
    params:
        options = "--max-depth 100 --min-BQ 15"
    conda:
        "envs/bcftools.yaml"
    wrapper:
        "v6.1.0/bio/bcftools/mpileup" 

rule bcftools_call:
    input:
        pileup = "results/variants/{sample}.bcf"
    output:
        calls = "results/variants/{sample}.vcf.gz"
    log:
        "logs/bcftools/call/{sample}.log"
    params:
        uncompressed_bcf = False,
        caller = "-m",
        extra = "--ploidy 1 --prior 0.001"
    conda:
        "envs/bcftools.yaml"
    wrapper:
        "v7.1.0/bio/bcftools/call"  
