## GISA-seq


#### Description
A pipeline for NGS analysis of Genome-wide Integration Site of AAV-Donor identified by sequencing (GISA-seq).

#### Required tools&packages
trimmomatic >= 0.39
flash >= 1.2.11
bwa >= 0.7.17-r1188
samtools >= 1.20
bedtools >= v2.30.0
perl >= v5.30.3 
python >= 3.8
annotatePeaks (http://homer.ucsd.edu/homer/ngs/annotation.html)

```shell
usage: GISAseq_pipeline.py [-h] -g GENOME -hr HOMOARM -s SITE [-ap ADAPTORS] [-tp TRIM_PARS] [-a ANNOTATION] [-gn GENOME_NAME] [-pl PLPROGRAM]

optional arguments:

  -h, --help            show this help message and exit
  -r2 READ2, --read2 READ2
                        [Required] Read 2 sequence file, .fastq or .fastq.gz
  -g GENOME, --genome GENOME
                        [Required] Indexed artificial genome file, .fasta
  -bc BARCODE, --barcode BARCODE
                        [Required] Downstream 8bp barcode after FWD primer.
  -hr HOMOARM, --homoarm HOMOARM
                        [Required] 3'HR between AAV_donor and target region.
  -s SITE, --site SITE  [Required] On-target integration site position. The format is chr:number. (e.g. chr9:46230897)

  -ap ADAPTORS, --adaptors ADAPTORS
                        [Optional] A FASTA file containing all adaptors need to be trimmed.
  -tp TRIM_PARS, --trim_pars TRIM_PARS
                        [Optional] Arguments and parameters for trimmomatic.
  -a ANNOTATION, --annotation ANNOTATION
                        [Annotation] Annoteated file for input genome, .GTF format.
  -gn GENOME_NAME, --genome_name GENOME_NAME
                        [Annotation] Name of the input genome, default is 'mm10'
  -pl PLPROGRAM, --plprogram PLPROGRAM
                        [Annotation] where is the annotatePeaks.pl

