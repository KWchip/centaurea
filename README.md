Hi there, this is a brief documentation on how I processed the RADseq data into `.vcf` files. 
# File structure
When using `snakemake` it is important to have an organised project `root` folder with their dedicated use case, below is an example:
## At a glance
```bash
.
├── config
├── dag.svg
├── data
├── envs
├── logs
├── README.md
├── reference
├── results
├── scripts
├── Snakefile
└── tmp
```
- `config`: contain configuration file like `config.yaml` and `config_small.yaml`, meaning the name of the input file and reference genome files are in a `yaml` format. such as:
`config_small.yaml`
```yaml
samples:
 - Spain20407
 - Australia60513
 - Chile40307

reference: "reference/centaurea_genomic.fna"
```
- `dag.svg`: This is a flowchart generated after we use `--dry-run` to test our `snakemake` (ignore for now, will only be important later on).
- `data`: This is where all the `{sample}_x_fq.gz` located. Basically the raw read files. 
- `envs`: This is where all the `bioconda` environment dependencies file lies (e.g. `bwa.yaml`, `samtools.yaml`, etc), usually, we will download each tools in a separate `conda` environments so that we can minimise conflicting dependencies. Example `conda` environment file:
`bwa.yaml`
```yaml
name: bwa
channels:
  - bioconda
  - conda-forge
dependencies:
  - _libgcc_mutex=0.1=conda_forge
  - _openmp_mutex=4.5=2_gnu
  - bwa=0.7.19=h577a1d6_1
  - bzip2=1.0.8=h4bc722e_7
prefix: /home/dell/Bioinformatics/Miniforge/envs/bwa
```
- `logs`: This is where you store all the log files, mostly useful when you encounter errors and are going on a debugging frenzy. It is recommended but not a must. 
- `README.md`: This is where you document things that users should know about this project root folder when they see one. Something like a welcoming poster, or an instruction on how to understand your project folder. 
- `reference`: Here is where we store our `.fna` reference genome file (in the beginning of the project) and the mapped files like `.pai`, `.fai`, etc. 
- `results`: This is where we store our final results. YIPEEE!!!
- `scripts`:This is where we store our utility coding scripts if there are any. (Something like housekeeping codes to make config file etc, not necessarily needed)
- `Snakefile`: The final boss
- `tmp`: This is where you keep all the experimental folder and files that you think you may need them later, but if deleted, it's all good nonetheless. 
## In-depth
So, when you are familiar with the above, you can now have a glance on this sub directory structure and get a gist of how the files are organised. 
```bash
config
├── config_small.yaml
└── config.yaml
envs
├── bcftools.yaml
├── bwa.yaml
├── fastqc.yaml
├── sambamba.yaml
├── samtools.yaml
└── trim_galore.yaml
reference
├── centaurea_genomic.fna
├── centaurea_genomic.fna.amb
├── centaurea_genomic.fna.ann
├── centaurea_genomic.fna.bwt
├── centaurea_genomic.fna.fai
├── centaurea_genomic.fna.pac
└── centaurea_genomic.fna.sa
scripts
├── config
│   └── config.yaml
└── make_config.sh
Snakefile  
tmp
README.md 
```
- having `config.yaml` and `config_small.yaml` is very handy because when you are writing your `Snakefile`, you may not want to repeat running the entire project data. Therefore, we can create a smaller subset of the data like a few samples to run through the entire pipeline to test for bug and reproducibility.
#  Setting up `conda` environment
Usually, i will create a `conda` environment for each specific tools (not a recommendation, just a user-to-user preference). So below is how I would create one:
```bash
conda create -n <env name> -c bioconda -c conda-forge -c defaults <env name> --yes
```
For example, if I want to download `bwa`
```bash
conda create -n bwa -c bioconda -c conda-forge -c defaults bwa --yes
```
After that, I will create a `bwa.yaml` on `envs/` with the code below:
```bash
conda env export --name <env_name> > <env>.yaml
```
For example,
```bash
cd envs
conda env export --name bwa > bwa.yaml
```
You can check how they look like by running this command,
```bash
cat bwa.yaml
```
You might see something like this, if not, something is wrong!!! (Oh no)
```
name: bwa
channels:
  - bioconda
  - conda-forge
dependencies:
  - _libgcc_mutex=0.1=conda_forge
  - _openmp_mutex=4.5=2_gnu
  - bwa=0.7.19=h577a1d6_1
  - bzip2=1.0.8=h4bc722e_7
```
# Snakefile
This is the trickiest part of the entire project. 
Please note that the `Snakefile` I share here may not be the best/optimal, it just works. 
```python
configfile: "config/config_small.yaml"
SAMPLES = config["samples"]
REF     = config["reference"]

rule all:
    input:
        expand("results/variants/{sample}.vcf.gz", sample=SAMPLES),
        expand("results/variants/{sample}.vcf.gz.tbi", sample=SAMPLES),


#########################################################
# 2. ADAPTER/TRIM
#########################################################
rule trim_galore_pe:
    input:
        ["data/{sample}.1.fq.gz", "data/{sample}.2.fq.gz"]
    output:
        fasta_fwd = "results/trimmed/{sample}_val_1.fq.gz",
        report_fwd = "results/trimmed/reports/{sample}_R1_trimming_report.txt",
        fasta_rev = "results/trimmed/{sample}_val_2.fq.gz",
        report_rev = "results/trimmed/{sample}_R2_trimming_report.txt"
    log:
        "logs/trim_galore/{sample}.log"
    threads: 4
    params:
        extra = "--illumina -q 20"
    conda:
        "envs/trim_galore.yaml"
    wrapper:
        "v7.1.0/bio/trim_galore/pe"


#########################################################
# 3. BUILD REFERENCE INDEX
#########################################################
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
        samtools faidx {input.ref} 2> {log}
        """

#########################################################
# 4. ALIGNMENT
#########################################################
rule bwa_mem:
    input:
        reads = ["results/trimmed/{sample}_val_1.fq.gz",
                 "results/trimmed/{sample}_val_2.fq.gz"],
        ref = REF,
        idx = multiext(REF, ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        "results/aligned/{sample}.bam"
    log:
        "logs/bwa/mem/{sample}.log"
    threads: 8
    params:
        # ensure the RG line is expanded per-sample
        extra = lambda wildcards: f"-R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA'",
        sorting = "samtools",
        sort_order = "coordinate"
    conda:
        "envs/bwa.yaml"
    wrapper:
        "v1.0.0/bio/bwa/mem"

#########################################################
# 4b. MARK DUPLICATES (important for eukaryotic genomes)
#########################################################
rule mark_duplicates_samtools:
    input:
        bam = "results/aligned/{sample}.bam"
    output:
        dedup_bam = "results/aligned/{sample}.dedup.bam",
        dedup_bai = "results/aligned/{sample}.dedup.bam.bai"
    log:
        "logs/samtools/markdup/{sample}.log"
    threads: 4
    conda:
        "envs/samtools.yaml"
    shell:
        r"""
        set -o pipefail
        # name-sort -> fixmate -> coord-sort -> markdup -> index
        samtools sort -n -@ {threads} {input.bam} -T results/aligned/{wildcards.sample}.ns -o results/aligned/{wildcards.sample}.name_sorted.bam 2>> {log}
        samtools fixmate -m results/aligned/{wildcards.sample}.name_sorted.bam results/aligned/{wildcards.sample}.fixmate.bam 2>> {log}
        samtools sort -@ {threads} results/aligned/{wildcards.sample}.fixmate.bam -T results/aligned/{wildcards.sample}.cs -o results/aligned/{wildcards.sample}.pos_sorted.bam 2>> {log}
        samtools markdup -r results/aligned/{wildcards.sample}.pos_sorted.bam {output.dedup_bam} 2>> {log}
        samtools index -@ {threads} {output.dedup_bam} 2>> {log}
        # optional cleanup of intermediates
        rm -f results/aligned/{wildcards.sample}.name_sorted.bam results/aligned/{wildcards.sample}.fixmate.bam results/aligned/{wildcards.sample}.pos_sorted.bam
        """

#########################################################
# 5a. MEAN DEPTH (empirical measurement after deduplication)
#########################################################
rule samtools_mean_depth:
    input:
        bam = "results/aligned/{sample}.dedup.bam"
    output:
        mean_txt = "results/coverage/{sample}.mean_depth.txt"
    log:
        "logs/coverage/{sample}.log"
    threads: 2
    conda:
        "envs/samtools.yaml"
    shell:
        r"""
        # compute genome-wide mean depth (includes positions with depth 0)
        samtools depth -a {input.bam} | awk '{{sum+=$3; cnt+=1}} END {{ if (cnt>0) printf("%.6f\n", sum/cnt); else print 0 }}' > {output.mean_txt} 2>> {log}
        """

#########################################################
# 5b. VARIANT CALLING (bcftools mpileup with dynamic --max-depth)
#########################################################
rule bcftools_mpileup:
    input:
        ref   = REF,
        faidx = REF + ".fai",
        bam   = "results/aligned/{sample}.dedup.bam",
        bai   = "results/aligned/{sample}.dedup.bam.bai",
        mean  = "results/coverage/{sample}.mean_depth.txt"
    output:
        pileup = "results/variants/{sample}.bcf"
    log:
        "logs/bcftools/mpileup/{sample}.log"
    threads: 4
    conda:
        "envs/bcftools.yaml"
    shell:
        r"""
        set -o pipefail
        mean=$(cat {input.mean} | tr -d '[:space:]')
        # guard against empty/zero mean (failed sample)
        if [ -z "$mean" ] || awk "BEGIN{{exit($mean<=0)}}"; then
            mean=0
        fi
        # compute max-depth = max(100, int(mean * 5))
        maxd=$(awk -v m="$mean" 'BEGIN {{ md = int(m*5); if (md < 100) md = 100; print md }}')
        echo "Sample {wildcards.sample} mean_depth=$mean, using --max-depth=$maxd" >> {log}
        # use mapping quality filter (-q) and base quality (-Q). Output compressed BCF.
        bcftools mpileup -f {input.ref} -q 20 -Q 15 --max-depth $maxd -Ob -o {output.pileup} {input.bam} 2>> {log}
        """

rule bcftools_call:
    input:
        pileup = "results/variants/{sample}.bcf"
    output:
        calls = "results/variants/{sample}.vcf.gz"
    log:
        "logs/bcftools/call/{sample}.log"
    params:
        caller = "-m",
        extra = "--ploidy 2 --prior 0.001"
    threads: 2
    conda:
        "envs/bcftools.yaml"
    wrapper:
        "v7.1.0/bio/bcftools/call"

rule bcftools_index:
    input:
        "results/variants/{sample}.vcf.gz"
    output:
        "results/variants/{sample}.vcf.gz.tbi"
    log:
        "logs/bcftools/index/{sample}.log"
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools index --tbi {input} 2> {log}"  
```
I recommend using `snakemake` wrapper as much as possible whenever available because they provide optimal reproducibility across platforms. 

Thank you,
Best, 
~Chip Kam Weng
