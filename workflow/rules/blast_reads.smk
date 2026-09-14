rule blast_reads:
    input:
        merged_prosyn_reads= scratch_dict["prosyn_reads"]["extracted_reads"] / "{sample}_merged.fasta", 
        diamond_db = Path(config["input"]["diamond_file"]),
    output:
        diamond_out = scratch_dict["diamond_blast"] / "{sample}.tsv",
    conda:
        "../envs/diamond-blast.yaml"
    shell:
        """
        if [ -s {input.merged_prosyn_reads} ]; then
            diamond blastx \
                --query {input.merged_prosyn_reads} \
                --db {input.diamond_db} \
                --out {output.diamond_out} \
                --threads {resources.cpus_per_task} \
                --outfmt 6 qseqid sseqid pident nident length qstart qend sstart send evalue bitscore \
                --max-target-seqs 1
        else
            touch {output.diamond_out}
        fi
        """
