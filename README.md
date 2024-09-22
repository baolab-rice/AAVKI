## GISA-seq


#### Description
A pipeline for NGS analysis of Genome-wide Integration Site of AAV-Donor identified by sequencing (GISA-seq).


```shell
usage: AAVKI_pipeline.py [-h] -g GENOME -hr HOMOARM -s SITE [-ap ADAPTORS] [-tp TRIM_PARS] [-v]

AAVKI pipeline is for quantification of CRIPSR editing outcomes with AAV integration. This pipeline
needs artificial genome with AAV vector used in experiment. For the artificial genome preparation
details, see the part 'Artifical Genome Prep' in our github.

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
  -a ANNOTATION, --annotation ANNOTATION
                        [Required] Annoteated file for input genome, .GTF format.
  -ap ADAPTORS, --adaptors ADAPTORS
                        [Optional] A FASTA file containing all adaptors need to be trimmed.
  -tp TRIM_PARS, --trim_pars TRIM_PARS
                        [Optional] Arguments and parameters for trimmomatic.
