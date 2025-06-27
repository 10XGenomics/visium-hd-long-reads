# Enabling Visium HD 3' + Long-Reads Analysis with Space Ranger

Amplified cDNA from the 10x Genomics Visium HD 3’ assay can be adapted for long-read sequencing technologies such as Pacific Biosciences (PacBio) and Oxford Nanopore Technologies (ONT). To process Visium HD 3’ long-read data, additional steps are required to prepare long-read sequences for the Space Ranger pipeline and to assign the Visium HD barcodes (BCs), unique molecular identifiers (UMIs), and spatial coordinates back to the original long-reads.

This repo walks through how to assign UMIs and corrected barcodes to Visium HD 3’ long-read sequencing data using custom python scripts and Space Ranger (version 4.0 and above). The overall approach is to first generate synthetic short paired-end reads from the original long reads using a custom pre-processing script [long_reads_to_10x_paired.py](https://github.com/10XGenomics/visium-hd-long-reads/blob/main/long_reads_to_10x_paired.py). These synthetic short reads can be used as input to Space Ranger `count`. A second custom script [add_10x_bam_tags.py](https://github.com/10XGenomics/visium-hd-long-reads/blob/main/add_10x_bam_tags.py) then transfers the Visium HD UMIs and corrected barcodes back to the long-read data to be followed by downstream analysis with the long-read sequencing provider software tools.

## Input Data

### Oxford Nanopore Technologies (ONT)

If using ONT, they have provided a conda package called [`percula`](https://anaconda.org/nanoporetech/percula) to preprocess ONT long-read BAMs for Space Ranger that will also be compatible to re-enter their [`wf-single-cell`](https://github.com/epi2me-labs/wf-single-cell) workflow for downstream analysis. See their support documentation and usage of `percula` [here](https://epi2me.nanoporetech.com/epi2me-docs/tools/percula/).

### Pacific Biosciences (PacBio)

Provide the Kinnex segmented BAM ([S-reads](https://skera.how/read-segments.html)) output from SMRT analysis software, or following `skera split` if using command-line workflow. More information on using `skera` can be found [here](https://skera.how/).

## Requirements
```
python>=3.11.*
python-edlib==1.3.*
pysam==0.22.*
```
Dependencies can be installed in a conda environment as follows:
```
$ conda create -n <env_name> -f requirements.txt -c conda-forge
$ conda activate <env_name>
```

# Step 1: Pre-Processing PacBio Data for Analysis with Space Ranger
Use the [`long_reads_to_10x_paired.py`](https://github.com/10XGenomics/visium-hd-long-reads/blob/main/long_reads_to_10x_paired.py) python pre-processing script to generate synthetic short paired-end reads from your long-read BAM file.
Example:
```
$ python3 long_reads_to_10x_paired.py \
	--bam <input_long_reads_bam> \
	--sample_name <sample_name> \
	[--optional params]
```
__Minimal Required Inputs__  
* `–-bam`: The long_reads.bam from Visium HD 3’ cDNA library, or a long-reads FASTQ (or `--fastq` <.fastq file>. Note that while FASTQ files can be used for the pre-processing step, we recommend converting to BAM because it is required for the post-processing step.)
* `--sample_name`: The desired prefix for all output files.

__Outputs__ 
* `*S1_L001_R1_001.fastq.gz, *S1_L001_R2_001.fastq.gz` (paired-end R1/R2 reads in FASTQ format)
* `summary.tsv.gz` (adapter configuration and location within each read)
* `configs.json` (count summary of adapter sub-read configurations)

_Optional parameters:_
* `--fastq`: (optional alternative to using the BAM) The long read FASTQ file to be used as input (note that post-processing script requires the raw long read data in BAM format, see instructions above for converting FASTQ to BAM format)
* `--threads`: The number of threads to use, default = `12`
* `--compress`: Compress the R1/R2 FASTQ outputs, default = `true`
* `--chunk_size`: The chunk size for multiprocessing, default = `2500`
* `--r1_size`: The number of bases after the adapter to include in R1 FASTQ file, default = `43`
* `--r2_size`: The number of bases after the adapter to include in R2 FASTQ file, default = `200`. Minimum read length recommendation is 75, maximum is 250.
* `--min_id`: The minimum sequence identity score to call an adapter, default = `0.8`

# Step2: Space Ranger Analysis
Use the R1/R2 paired-end FASTQs from the pre-processing script as input for `spaceranger count`, see additional documentation on our support site for processing Visium HD ‘3 data [here](). 

__It’s important to set the `--create-bam=true` parameter in your command__, as the `possorted_genome_bam.bam` produced by `spaceranger count` will be required for the post-processing step to map the corrected barcode and UMI sequences back to the original long-read BAM.

Example:
```
$ spaceranger count --id="HD_Adult_Mouse_Brain" \
      --transcriptome=refdata-gex-mm10-2024-A \
      --fastqs=datasets/HD_Adult_Mouse_Brain_fastqs \
      --image=datasets/HD_Adult_Mouse_Brain_image.tif \
      --slide=V19L01-041 \
      --area=A1 \
      --localcores=16 \
      --localmem=128 \
      --create-bam=true
```

# Step 3: Post-Processing Space Ranger Outputs
The post-processing script [`add_10x_bam_tags.py`](https://github.com/10XGenomics/visium-hd-long-reads/blob/main/add_10x_bam_tags.py) tags the raw, full-length long-read input BAM with the corrected Visium HD 3’ barcodes and UMIs. You can then use this file for additional long-read processing, alignment, and other secondary analysis workflows. The following standard BC/UMI BAM tags are used (see [10x BAM tag documentation](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/space-ranger-bam) for additional details): 
* `CB`: Corrected barcode (Note: reads that cannot be assigned a corrected barcode will not have a `CB` tag.)
* `CR`: Uncorrected barcode
* `UB`: Corrected UMI
* `UR`: Uncorrected UMI

Usage:
```
$ python3 add_10x_bam_tags.py \
	--sr_bam <possorted_genome_bam.bam> \
	--lr_bam <long-reads.bam> \
	[--optional params]
```
__Minimal Required Inputs__
* `--sr_bam`: `possorted_genome_bam.bam` from the Space Ranger outputs
* `--lr_bam`: Original long-read BAM used in pre-processing

__Outputs__
* `<lr_bam>.tagged.bam` (adds BC/UMI tags to long-reads BAM `CR`,`CB`,`UR`,`UB`)
* `*.spatial_barcodes.csv.gz` (CSV containing the original long-read name identifiers and the BC/UMI tags). Example:
   ```
   read_name,uncorrected_barcode,corrected_barcode,uncorrected_umi,corrected_umi
   m84039_250124_023230_s2/151066185/ccs/12245_15157,GTCTGCATCTGCCCTGCATTAATGCATCAG,s_002um_02077_01449-1,CTGGGACGA,CTGGGACGA
   m84039_250124_023230_s2/187635340/ccs/2043_2623,GCAGCTATGCAGGTAGTATCCACGGCATCG,s_002um_00964_02405-1,CAATGCATA,CAATGCATA
   m84039_250124_023230_s2/146084225/ccs/2498_3082,GCAGCTATGCAGGTAGTATCCACGGCATCG,s_002um_00964_02405-1,CAATGCATA,CAATGCATA
   ```
_Optional parameters:_ 
* `--chunk_size`: Chunk size for multiprocessing, default = `20000`
* `--out_dir`: The desired output path, default = current working directory 
* `--overwrite`: Option to overwrite the raw long-reads BAM instead of generating a new `*.tagged.bam`, default = `false`

# Continuing Your Visium HD 3’ Long Read Analysis
With long-read identifiers now linked to their BC and UMI tags, the `tissue_position.parquet` file can be used to assign spatial coordinates to the identifiers. This file directly maps the center of each barcode to x and y pixel coordinates (from the `pxl_col_in_fullres` and `pxl_row_in_fullres` columns) in the full-resolution microscope image. Subsequently, the `*.spatial_barcodes.csv.gz` file connects these barcodes (and by extension their pixel coordinates) to their respective long-read identifiers, which then enables the visualization of each identifier's spatial location.









