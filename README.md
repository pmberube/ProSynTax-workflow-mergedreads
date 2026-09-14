# ProSynTax Workflow: taxonomic classification of *Prochlorococcus* and *Synechococcus*


## Introduction
Here we present a workflow that accompanies ProSynTax - a curated protein sequence dataset aimed at enhancing the taxonomic resolution of *Prochlorococcus* and *Synechococcus* classification. ProSynTax includes proteins from 1,260 genomes of *Prochlorococcus* and *Synechococcus*, including single amplified genomes, high-quality draft genomes, and newly closed genomes. Additionally, ProSynTax incorporates proteins from 41,753 genomes of marine heterotrophic bacteria, archaea, and viruses to assess microbial and viral communities surrounding *Prochlorococcus* and *Synechococcus*. This resource enables accurate classification of picocyanobacterial clusters/clades/grades in metagenomic data – even when present at 0.15% of reads for *Prochlorococcus* or 0.03% of reads for *Synechococcus*. 

![Phylogenetic Tree](docs/images/figure1.svg "Phylogenetic Tree")

## Publication 
This workflow and the accompanying ProSynTax dataset is described in: 

         Coe, A., Mullet, J.I., Vo, N.N. et al. A curated protein dataset for taxonomic classification of Prochlorococcus and Synechococcus in metagenomes. Sci Data 12, 1895 (2025). https://doi.org/10.1038/s41597-025-06164-5

## About this fork

This is a modified fork of [jamesm224/ProSynTax-workflow](https://github.com/jamesm224/ProSynTax-workflow), adapted to classify **merged reads** instead of paired-end reads.

### Why merged reads

Kaiju and DIAMOND BLASTX each classify one sequence at a time. In the upstream workflow both mates are passed to Kaiju, but only the forward read of each pair is carried through to the CyCOG alignment step, so the reverse read contributes nothing to the genome equivalent calculation.

Merging overlapping read pairs before classification gives a single query per fragment that is both longer and more accurate. The overlapping region is sequenced twice, so merging resolves discordant base calls and raises the effective quality of that region, and the merged fragment is longer than either mate alone. Longer, higher-quality queries yield longer and more reliable alignments against ProSynTax and CyCOG6, and alignment length is exactly what the genome equivalent normalization is computed from.

Reads must therefore be quality-filtered, adapter-trimmed and merged **before** this workflow is run. The workflow takes one merged FASTQ per sample as its input.

### What differs from upstream

| Area | Upstream | This fork |
|---|---|---|
| Input | Paired-end FASTQ; `forward read` + `reverse read` columns | One merged FASTQ; `merged_read` column |
| Trimming | BBDuk adapter + quality trimming, `minlen=25` | BBDuk length filter only, `minlen=150` |
| Kaiju | Both mates passed with `-i r1 r2` | Single merged file |
| DIAMOND query | Forward read only (`_fwd.fasta`) | Merged read (`_merged.fasta`) |
| Samples with no Pro/Syn reads | Run aborts | Empty result written, run continues |
| Normalization denominator | Hardcoded constant | Summed from `inputs/cycog_len.tsv` at runtime |
| Snakemake | v7 with a `--cluster` submission string | v9.19 with the SLURM executor plugin |
| Job submission | `sbatch run_classify_smk.sbatch` | `snakemake --profile profile`; script removed |
| BBMap | 39.18 | 39.37 |

Adapter and quality trimming are dropped from the workflow because both are assumed to have happened before merging. If that is not true of your data, restore those BBDuk parameters in `workflow/rules/trim_reads.smk`.

### Known gaps in this fork

- `inputs/example/` still contains paired-end example data, and `inputs/example/make_samples_tsv.py` still writes the upstream three-column schema, so the bundled example does not run against this version.
- The workflow diagram in `docs/images/figure2.svg` depicts the upstream paired-end pipeline.

## Bug Fix June 2026. 
A bug was identified that resulted in undercounting genome equivalents by a factor of 3. In the normalization by average cycog length, the denominator was expressed in terms of nucleotides rather than amino acids (script "normalize_all_cycog.py"). This calculation is updated in the bug fix branch "GE_bugFix". 

This fork carries that fix, and goes one step further: the denominator is not hardcoded at all. `normalize_all_cycog.py` sums the mean lengths in `inputs/cycog_len.tsv` at runtime, so the denominator and the CyCOG filter list are always drawn from the same rows of the same file and cannot drift apart.

## Restored 424th core CyCOG

`inputs/cycog_len.tsv` as distributed upstream held only **423** of the 424 single-copy core CyCOGs. The missing entry is **CyCOG_60001271** (Protein of unknown function, DUF2518; 387 genes; mean length 157.73385012919897 aa). The same record is damaged in the Zenodo `average_cycog_length.csv`, whose final line reads `CyCOG_60001271,` with no value and no trailing newline — both files were truncated at the same last record.

The core set was reconstructed to confirm this. The 423 published CyCOGs are present in exactly one copy in exactly 93 genomes (67 *Prochlorococcus*, 26 *Synechococcus*); asking in turn which CyCOGs are single-copy in all 93 of those genomes returns exactly 424, being the published 423 plus CyCOG_60001271, and matching the 424 identifiers listed in the Zenodo CSV. Mean lengths were recomputed from `CyCOG6_clean_prokkafmt_simple.faa` over all genes assigned to each CyCOG; all 423 published values were reproduced exactly (maximum absolute difference 0.0).

This fork ships the complete 424-entry table, which raises the denominator from 119801.5464 to 119959.2802. **Restoring the entry does not systematically shift genome equivalents.** `normalize_all_cycog.py` filters reads to exactly the CyCOGs listed in this file and then divides by the summed mean length of those same CyCOGs, so the estimator is self-consistent at any marker-set size and a 423-marker table is no more biased than a 424-marker one. Adding the marker raises the denominator by 157.73 aa, and for a clade present at G genome equivalents it simultaneously raises the numerator by roughly G x 157.73 aa as reads best-matching CyCOG_60001271 now pass the filter. The two changes cancel.

CyCOG_60001271 is a single-copy core gene, present in every genome of the reference set, so it recruits reads wherever *Prochlorococcus* or *Synechococcus* are present in a sample. The cancellation therefore holds in practice as well as in principle, and per-sample values move only slightly, in either direction, with read sampling on one additional marker. The complete table is shipped because it is correct and because one more marker contributes marginally more sequence to each estimate, not because it changes results.

## Table of Contents
* [About this fork](#about-this-fork)
  * [Why merged reads](#why-merged-reads)
  * [What differs from upstream](#what-differs-from-upstream)
  * [Known gaps in this fork](#known-gaps-in-this-fork)
* Setting up the Workflow
  * [Installing the ProSynTax Workflow](#installing-the-prosyntax-workflow) 
  * [Installing ProSynTax](#installing-prosyntax)
  * [Installing Dependencies](#installing-dependencies)
  * [Edit Workflow Parameters](#edit-workflow-parameters)
* Running the Workflow
  * [Submitting Main Script](#submitting-main-script)
  * [Troubleshooting Guide](#troubleshooting-guide)
* [Pipeline Workflow](#pipeline-workflow)
* Results
  * [Raw Counts Output](#1-raw-counts-output)
  * [Normalized Genome Equivalent Output](#2-normalized-genome-equivalent-output)
* [Limit of Detection Filtering](#limit-of-detection-filtering)
* [Intermediate Files](#intermediate-files)
* [How to build your own version of ProSynTax](#how-to-build-your-own-version-of-prosyntax)

## Setting up the Workflow
### Installing the ProSynTax Workflow
1. Clone the ProSynTax-workflow Github repository into your working directory:  

       # (optional) create a new project directory
       mkdir my_classification_project
       cd my_classification_project  # change directory into project 

- Create a copy of the workflow into your current path:  

       git clone https://github.com/pmberube/ProSynTax-workflow-mergedreads/

### Installing ProSynTax
All files associated with this workflow can be downloaded from the [Zenodo repository](https://doi.org/10.5281/zenodo.14889680) (DOI 10.5281/zenodo.14889680). Make sure to use the most up to date version of Zenodo for downloading files (v2 as of July 2025).

Download the following **required** files into a directory on your machine: 
1. ProSynTax_nodes.dmp
2. ProSynTax_names.dmp
3. ProSynTax_v1.1.fmi.bz2
4. CyCOG6.dmnd

- **Notes**:
  - Make sure to unzip ProSynTax_v1.1.fmi.bz2 prior to usage `bzip2 -dk ProSynTax_v1.1.fmi.bz2`
  - The path to these files will be needed later in the `inputs/config.yaml` file in step 1 of section [Edit Workflow Parameters](#edit-workflow-parameters). 

### Installing Dependencies
1. Install Mamba following instructions on the [official Mamba documentation](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). Or follow the steps below for installation into Linux system: 

- Download miniforge (contains mamba) from Github (for Linux):  

       wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

- Run the installer and set up mamba on your computer. :  

       # when prompted, press ENTER or enter YES
       bash Miniforge3-Linux-x86_64.sh 

- Once installation is complete, open a new terminal window  

       # check that installation worked
       mamba -V  # prints version 

- (Optional) This makes it so the base conda environment won’t automatically activate whenever you open a new terminal:  

       conda config --set auto_activate_base false 

- (Optional) Delete the installer  

       rm Miniforge3-Linux-x86_64.sh

2. Then, install [Snakemake v9.19](https://snakemake.readthedocs.io/en/stable/) and the SLURM executor plugin into a new mamba environment: 

       # create a new mamba environment named snakemake
       mamba create -c bioconda -c conda-forge -n snakemake snakemake=9.19 snakemake-executor-plugin-slurm
   
- (Optional) test that your installation worked: 

       mamba activate snakemake  # activate the environment 
       snakemake --version  # print snakemake version 
       mamba deactivate  # deactivate snakemake environment 

- **Note**: upstream uses Snakemake v7.32.4. The `profile/config.yaml` in this fork uses the executor plugin interface introduced in Snakemake 8 and will not work under v7. 
- **Why 9.19 is pinned**: this is the version this fork was developed and run against, and the `runtime` values in `profile/config.yaml` are matched to how 9.19 interprets them. Pinning keeps the workflow reproducible as run.
- **If you use a newer Snakemake instead**: the meaning of `runtime` changed in Snakemake 9.20, so the values in `profile/config.yaml` will request far more wall time than intended and your jobs may be rejected outright. Read the [warning about `runtime`](#a-warning-about-runtime) and adjust those values before your first run. Nothing else in the workflow is known to depend on the Snakemake version.

### Edit Workflow Parameters
#### 1. Edit ```inputs/config.yaml```:  
This is the experimental configuration file, where you will specify ProSynTax file paths, directory paths, and more. 

- **Required** edits:  
  - `nodes_file`: path to your installation of ProSynTax_nodes.dmp file for Kaiju (completed earlier in step [Installing ProSynTax](#installing-prosyntax)). 
  - `names_file`: path to your installation of ProSynTax_names.dmp file for Kaiju. 
  - `fmi_file`: path to your installation of ProSynTax_file.fmi file for Kaiju. 
  - `diamond_file`: path to your installation of CyCOG6.dmnd dataset file for BLAST. 
  - `genus_list`: list of genus you would like to extract raw read count for. 
    - Default: ['Synechococcus', 'Prochlorococcus', 'unclassified']
    - "unclassified": include reads that are labeled "unclassified" or "cannot be assigned to a (non-viral) genus" by Kaiju. 
    - Reads classified as other genus not listed will be summed into 1 group called "other_genus". 
    - Refer to [Results](#results) for more information. 
  - `scratch directory`: path to folder for storing intermediate files. 
    - Such as: trimmed read files, Kaiju outputs, and Blast outputs. 
    - Shipped as a `/path/to/` placeholder; set it to an absolute path on scratch storage. 
  - `results directory`: path to directory for storing final output files: "summary_read_count.tsv" and "normalized_counts.tsv".
    - Also a `/path/to/` placeholder that must be set.

#### 2. Create ```inputs/samples.tsv```: 
This file should contain metadata for your samples, one row per sample. 

- **Required** columns: 
  - `sample`: unique name for sample. 
  - `merged_read`: path to the merged read FASTQ file for that sample. 
- Optional: 
  - Feel free to add any other sample metadata columns as they will not impact the workflow. 

Tab-delimited, for example:

    sample	merged_read
    SAMPLENAME	/PATH/TO/MERGED/FASTQ

- **Note**: the helper script [inputs/example/make_samples_tsv.py](inputs/example/make_samples_tsv.py) has **not** been updated for this fork and still writes the upstream `forward read` / `reverse read` columns. Write `samples.tsv` yourself, or adapt that script.

#### 3. Edit ```profile/config.yaml```:  
This file contains Snakemake specifications and compute resource (HPC) specifications. For tips on determining your HPC cluster resource specification for this section, visit [HPC Resource Tips](docs/readme_extras/resource_tips.md).

This fork drives SLURM through the Snakemake executor plugin rather than the `--cluster` submission string used upstream, so the resource keys differ from those described in the upstream README. The shipped values are specific to the MIT Engaging cluster and must be changed.

- **Required** edits:
  - `jobs`: number of jobs you would like to run at once on the compute cluster. 
  - `default-resources: slurm_partition`: partition name to submit jobs to. 
    - To check what partitions you have access to: `sinfo`.
  - `default-resources: slurm_account`: account to charge the jobs to.
  - `default-resources: runtime`: wall time per job. **Give this a unit** — see the warning below.
  - `default-resources: mem_mb`: memory in MB for rules not listed under `set-resources`. 
  - `set-resources: <rule>: cpus_per_task`: number of cores to allocate for multi-threaded rules. 
  - `set-resources: <rule>: mem_mb`: memory in MB for individual rules. 
- Optional edits:
  - `default-resources: slurm_extra`: additional flags passed through to sbatch. The shipped value is `'--requeue'`, which suits a preemptable partition.
  - `conda-prefix`: path to a previous conda installation, located in the hidden `.snakemake` directory.
    - If you have run this pipeline before, you can save time on newer runs by referencing previous conda installations.  
    - For example, if you ran this workflow in `/home/my_username/project1/` and would like to run the workflow on another project, the "conda-prefix" path for the newer project would be: "/home/my_username/project1/.snakemake/conda"

##### A warning about `runtime`

How Snakemake reads a bare number in `runtime` **changed in Snakemake 9.20**, so the same profile means different things on either side of that release:

| Snakemake | `runtime=43200` is read as | Result |
|---|---|---|
| 9.19 and earlier | 43200 **seconds** | 12 hours |
| 9.20 and later | 43200 **minutes** | 30 days |

In 9.19 and earlier, a bare number was passed to `humanfriendly.parse_timespan`, which treats a unitless value as seconds, and the result was divided by 60 to give minutes. From 9.20 a decimal string short-circuits that path and is taken as minutes directly.

This is not hypothetical. Under Snakemake 9.19 with `runtime=720`, SLURM allocated **12 minutes** rather than 12 hours (720 / 60), killing `kaiju_run` partway through. The shipped value of 43200 is the correction for that behaviour: under 9.19 it gives the intended 12 hours.

**On any current Snakemake, the shipped value of 43200 requests 30 days per job and will be rejected by a partition with a shorter `MaxTime`** (`mit_preemptable` allows 2 days).

The install instructions pin Snakemake to 9.19 so that the shipped values mean what they say. If you deliberately move to a newer release, either recompute every `runtime` value in minutes, or give each one a unit — a string with units is parsed identically on every version, so this form is safe either way:

    - runtime="12h"

Check what your partition actually allows with `scontrol show partition <name>`.


## Running the Workflow
### Submitting Main Script
The executor plugin submits every rule as its own SLURM job, so this fork has no wrapper submission script; `run_classify_smk.sbatch` has been removed. Activate the environment and run Snakemake directly:

    mamba activate snakemake
    snakemake --profile profile

- Run `snakemake --profile profile -n` first for a dry run listing the jobs that would be submitted.
- Snakemake itself must stay alive while its jobs run, so start it under `tmux`, `screen`, or `nohup` for long runs.


Track Workflow Progress: 
- Snakemake reports the progress of each step (called a rule) on its own standard output. Redirect it to a file, or run under `tmux`/`screen`, if you want to read it back later. 
- The SLURM executor plugin writes per-job SLURM log files relative to the working directory; change this with `--slurm-logdir`. Logs for **successful** jobs are deleted on completion by default and failed jobs keep theirs, so pass `--slurm-keep-successful-logs` if you want all of them. 
- Upstream's `logs/{rule}/{wildcards}.err` layout came from the `--cluster` submission template and no longer applies. 
- Once all steps in the pipeline have completed, Snakemake's output ends with the following:
  ```
  N of N steps (100%) done
  Complete log: .snakemake/log/YYYY-MM-DDTXXXXXX.XXXXXX.snakemake.log
  ```

### Troubleshooting Guide
We have some tips for troubleshooting and debugging for common errors located [here](docs/readme_extras/debugging_tips.md). 


## Pipeline Workflow
Below is a description of the steps in ProSynTax workflow. 

![Workflow Overview](docs/images/figure2.svg "Pipeline Workflow")

**1. Quality Control:**  
  - This workflow expects reads that have already been quality-filtered, adapter-trimmed and merged, so the only remaining step is a length filter applied with [BBDuk](https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/) from the BBTools suite: 
    - `minlen=150`: discard merged reads shorter than 150bp 
  - Upstream instead performs adapter and quality trimming of raw paired-end reads at this step (`qtrim=rl`, `trimq=10`, `ktrim=r`, `k=23`, `mink=11`, `hdist=1`, `minlen=25`). Those parameters are removed here because merging is assumed to follow trimming. Restore them in `workflow/rules/trim_reads.smk` if that does not hold for your data.

**2. Taxonomic Classification:**  
  - Merged reads are classified against the ProSynTaxSB using [Kaiju v1.10.1](https://github.com/bioinformatics-centre/kaiju) with parameters: 
    - `-m 11`: Minimum match length
    - `-s 65`: Minimum match score in Greedy mode
    - `-E 0.05`: Minimum E-value in Greedy mode
    - `-x`: Enable SEG low complexity filter
    - `-e 5`: Number of mismatches allowed in Greedy mode

**3. Raw Read Count Summary:** 
  - Kaiju classification summaries are obtained using `kaiju2table` command. 
    - All output files are aggregated into final results table "summary_read_count.tsv" using `workflow/scripts/classification_summary.py`.   

**4. Clade Normalization:**  
  - Full taxon paths are added to Kaiju output in Step 2 using `kaiju-addTaxonNames` command with the parameters: 
    - `-p`: print the full taxon path instead of just the taxon name.
    - `-u`: omit unclassified reads (saves disk space compared to output of all reads).
  - FASTQ header name and full taxonomic classification of reads classified as "*Prochlorococcus*" and "*Synechococcus*" by Kaiju are extracted using the Bash command `grep`
     - Extracted FASTQ header names are used to pull the corresponding sequences out of the **merged** read file using [seqtk vr82](https://github.com/lh3/seqtk) for downstream BLAST. Upstream extracts from the forward read file only. 
  - Extracted FASTA sequences are compared against the CyCOGs listed in `inputs/cycog_len.tsv` using DIAMOND [Blastx v2.1.11](https://github.com/bbuchfink/diamond) sequence aligner with parameter: 
    - `max-target-seqs 1`: maximum number of target sequences per query to report alignments for. 
  - Using the custom Python script `normalize_all_cycog.py`, reads are normalized by clade following these steps: 
    - Filter for reads with hits to the CyCOGs listed in `inputs/cycog_len.tsv` 
    - Obtain sum of alignment length 
    - Divide that sum by the total CyCOG length, computed at runtime as the sum of column 2 of `inputs/cycog_len.tsv` (119959.2802 for the file as shipped) 
      - DIAMOND reports alignment length in amino acids and `cycog_len.tsv` is in amino acids, so numerator and denominator share units. See [Bug Fix June 2026](#bug-fix-june-2026) above. 
  - All normalized output files are aggregated into file results table `normalized_counts.tsv`. 
  - A sample in which Kaiju finds no *Prochlorococcus* or *Synechococcus* reads yields an empty entry rather than aborting the run.


## Results
This pipeline outputs 2 main output files (listed below), both located in "results directory" path specified in ```inputs/config.yaml``` file. 

### 1. Raw Counts Output
The output file `summary_read_count.tsv` contains the total classified raw read counts for the genus listed in ```inputs/config.yaml``` file.   
- Column Description:  
  - `sample_name`: name of sample
    - Same as the "sample" column in ```samples.tsv```.
  - `taxon_name`: name of genus 
    - Taxons in this column are defined by the genus list specified in `inputs/config.yaml` file. 
  - `reads`: number of reads Kaiju classified as specified genus in "taxon_name" 
  - `percent`: percent of this genus relative to all reads in the sample 

### 2. Normalized Genome Equivalent Output
The output file `normalized_counts.tsv` contains normalized genome equivalent for each *Prochlorococcus* and *Synechococcus* clade/subclade/ecotype in the sample.  
- Column Description:  
  - `sample_name`: name of sample; same as the "sample" column in ```samples.tsv```.
  - `genus`: genus classified by Kaiju ("*Prochlorococcus*" or "*Synechococcus*")
  - `clade`: ecotype/cluster/clade/grade classified by Kaiju 
  - `alignment_length`:  sum of alignment length of unique hits from Diamond Blast 
  - `genome_equivalents`: normalized abundance of classified ecotype/cluster/clade/grade in sample 
    - Refer to "Read Normalization" Step in [Pipeline Workflow](#pipeline-workflow) for more information on how this value was calculated 

## Limit of Detection Filtering
To assess the accuracy of *Prochlorococcus* and *Synechococcus* cluster/clade/grade classifications, we generated mock metagenome sequence datasets and conducted simulations following the methods outlined in Coe et al. (2025) to calculate misclassification rates. To maintain misclassification rates below either 10% or 5%, we determined the minimum *Prochlorococcus*:*Synechococcus* (or vice versa) ratios, as detailed below. For broader classifications at the cluster or ecotype level, a 10% misclassification rate is sufficient.  However, for more precise subcluster/clade/grade delineations, we recommend using a 5% false positive rate threshold.

Scripts used for genome subsetting and read simulation are located in the [validation](validation) directory. 

### *Prochlorococcus* Filtering Parameters
#### 5% Misclassification Parameters  
- *Prochlorococcus* Abundance > 0.15% 
  - Counts of *Prochlorococcus* reads out of *all* classified reads 
- Ratio of *Prochlorococcus*:*Synechococcus* > 0.40
  - Ratio calculated by counts of *Prochlorococcus* divided by counts of *Synechococcus*

#### 10% Misclassification Parameters
- *Prochlorococcus* Abundance > 0.08% 
  - Counts of *Prochlorococcus* reads out of all classified reads 
- Ratio of *Prochlorococcus*:*Synechococcus* > 0.23
  - Ratio calculated by counts of *Prochlorococcus* divided by counts of *Synechococcus*

### *Synechococcus* Filtering Parameters
#### 5% Misclassification Parameters  
- *Synechococcus* Abundance > 0.03% 
  - Counts of *Synechococcus* reads out of all classified reads 
- Ratio *Synechococcus*:*Prochlorococcus* > 0.20 
  - Ratio calculated by counts of *Synechococcus* divided by counts of *Prochlorococcus*

#### 10% Misclassification Parameters
- *Synechococcus* Abundance > 0.02% 
  - Counts of *Synechococcus* reads out of *all* classified reads 
- Ratio *Synechococcus*:*Prochlorococcus* > 0.10
  - Ratio calculated by counts of *Synechococcus* divided by counts of *Prochlorococcus*


## Visualization
For visualization of clade/grade/cluster composition, we have devised a list of HEX color codes, which are located in the [taxon_colors](docs/taxon_colors) directory. 

## Intermediate Files
This workflow produces several intermediate files that may be useful for additional analysis. Descriptions of these files can be found [here](docs/readme_extras/intermediate_files.md). 

## How to build your own version of ProSynTax
In the case that you wish to develop your own version of ProSynTax, we provided a smaller version of the dataset without any NCBI RefSeq Genomes (ProSynTax_v1.1_without_refseq.faa.bz2). This dataset contains the annotated *Prochlorcoccus*, *Synechococcus*, marine heterotrophs, and all of GORG-tropics. You can easily add additional genomes to this dataset. If you run into errors, ensure to update the names and nodes files with updated NCBI taxonomy.

To add new genomes:
1. Convert nucleotide sequences to protein sequences (ensure that no stop codons are present within the sequences).
2. For each genome ensure you create a unique identifier. They must be added to the end of each protein sequences.

            genome1-orf1_{identifier_number} -> i.e. >MIT0701-NODE1635_8000886
   
4. Update the names/nodes files to ensure the new genomes are included. It must include NCBI names/nodes files and the unique names/nodes files for ProSynTax genomes.
5. Add protein sequences from your new genomes to the ProSynTax_v1.1_without_refseq.faa.bz2 (on Zenodo).
6. Convert the new customized faa file into an fmi file following instructions at the Kaiju Github: [https://github.com/bioinformatics-centre/kaiju](https://github.com/bioinformatics-centre/kaiju#:~:text=to%20run%20Kaiju.-,Custom%20database,-It%20is%20also)

