rule kaiju_name_extract:
    input:
        kaiju_name = scratch_dict["classified_kaiju_read_output"] / "{sample}_names.out",
    output:
        read_name_file = scratch_dict["prosyn_reads"]["read_name"] / "{sample}.txt", 
        read_name_taxa_file = scratch_dict["prosyn_reads"]["read_name_classification"] / "{sample}_classification.txt", 
    shell:
        """
        # obtain name of reads whose classification contains Pro/Syn 
        grep -E "Prochlorococcus|Synechococcus" {input.kaiju_name} | cut -f2 > {output.read_name_file}

        # obtain name and full taxonomic classification of reads whose classification contains Pro/Syn 
        grep -E "Prochlorococcus|Synechococcus" {input.kaiju_name} | cut -f2,4 > {output.read_name_taxa_file}
        """

rule extract_fastq_reads:
    """
    Note: extracting from the merged reads fastq file. 
    """
    input: 
        r1 = scratch_dict["trimmed_reads"] / "{sample}_1_trimmed.fastq.gz",
        prosyn_read_name = scratch_dict["prosyn_reads"]["read_name"] / "{sample}.txt", 
    output:
        merged_prosyn_reads = scratch_dict["prosyn_reads"]["extracted_reads"] / "{sample}_merged.fastq", 
    conda:
        "../envs/seqtk.yaml"
    shell:
        """
        seqtk subseq {input.r1} {input.prosyn_read_name} > {output.merged_prosyn_reads}
        """

rule extracted_fastq_to_fasta: 
    """
    """
    input: scratch_dict["prosyn_reads"]["extracted_reads"] / "{sample}_merged.fastq", 
    output: scratch_dict["prosyn_reads"]["extracted_reads"] / "{sample}_merged.fasta", 
    conda: "../envs/seqtk.yaml"
    shell: "seqtk seq -A {input} > {output}"

