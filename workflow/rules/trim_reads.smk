rule run_trim_PE:
    input:
        r1 = lambda wildcards: SAMPLE_TABLE.loc[wildcards.sample, 'merged_read'],
    output:
        o1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq.gz",
    conda:
        "../envs/bbtools.yaml"
    shell:
        "bbduk.sh threads={resources.cpus_per_task} "
        "in={input.r1} "
        "out={output.o1} "
        "minlen=150 "
