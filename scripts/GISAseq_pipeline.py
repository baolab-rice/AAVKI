from GISAseq_QC_and_merge import run_QC_and_merge
from GISAseq_align import run_align, samtobed
from GISAseq_layered_filter import run_layered_filters
from GISAseq_OT_quantification import run_quantification, run_ON_quantification
from GISAseq_annotation import run_annotation
import re
import argparse

"""
Required programs: 
trimmomatic >= 0.39
flash >= 1.2.11
bwa >= 0.7.17-r1188
samtools >= 1.20
bedtools >= v2.30.0
perl >= v5.30.3 
python >= 3.8
"""

"""
This program should be run in the folder with sample.
"""

"""
Usage: python3 GISAseq_pipeline.py -r1 <read1.fq> -r2 <read2.fq> -g <artificial_genome.fa> -bc "NNNNNNNN" -ap <adaptors.fa>

For gene annotation, the -a, -gn and -pl flags are required,
"""


def main():
    # Quality control w/ raw reads and merge the reads
    # reversed read 1 & read 2?
    merged_reads = run_QC_and_merge(args.read2,args.adaptors,args.trim_pars,args.barcode)
    print("QC done.")

    # # Sequence alignment
    samfile = run_align(merged_reads,args.genome)
    print("Sequence alingment done.")

    # Layered filters
    OT_HDR_file,OT_NHEJ_file,ON_HDR_file,ON_NHEJ_file = run_layered_filters(samfile,int(args.homoarm),args.site)
    print("Filtering done.")

    # Quantification
    print("Quantification of OT_HDR")
    file_OT_HDR =  run_quantification(OT_HDR_file,args.site)
    print("Quantification of OT_NHEJ")
    file_OT_NHEJ = run_quantification(OT_NHEJ_file,args.site)
    print("Quantification of ON_HDR")
    run_ON_quantification(ON_HDR_file,args.site)
    print("Quantification of ON_NHEJ")
    run_ON_quantification(ON_NHEJ_file,args.site)
    print("Quantification done.")

    # Annotation
    if args.annotation != "":
        sam_OT_HDR = run_align(file_OT_HDR,args.genome)
        sam_OT_NHEJ = run_align(file_OT_NHEJ,args.genome)
    
        bed_OT_HDR = samtobed(sam_OT_HDR)
        bed_OT_NHEJ = samtobed(sam_OT_NHEJ)
    
        run_annotation(bed_OT_HDR,args.annotation, args.plprogram, args.genome_name)
        run_annotation(bed_OT_NHEJ,args.annotation,  args.plprogram, args.genome_name)
        print("Annotation done.")



if __name__ == '__main__':

    desc = "AAVKI pipeline is for quantification of CRIPSR editing outcomes \
        with AAV integration. This pipeline needs artificial genome with AAV vector \
        used in experiment. For the artificial genome preparation details, see the part \
        'Artifical Genome Prep' in our github."

    parser = argparse.ArgumentParser(description=desc)

    # Required arguments
    ## InputFile: Read1.fq.gz Read2.fq.gz artificial_genome.fa (indexed)
    # parser.add_argument('-r1','--read1',required=True,help="[Required] Read 1 sequence file, .fastq or .fastq.gz")
    parser.add_argument('-r2','--read2',required=True,help="[Required] Read 2 sequence file, .fastq or .fastq.gz")
    parser.add_argument('-g','--genome',required=True,help="[Required] Indexed artificial genome file, .fasta")
    
    ## Non-specific PCR filter
    parser.add_argument('-bc','--barcode',required=True,help="[Required] Downstream 8bp barcode after FWD primer.")

    ## The length of 3'HR
    parser.add_argument('-hr','--homoarm',required=True,help="[Required] 3'HR between AAV_donor and target region.")   

    ## Integrated site position
    parser.add_argument('-s','--site',required=True,help="[Required] On-target integration site position. The format is chr:number. (e.g. chr9:46230897)")

    # Optional arguments
    ## Trim in QC
    ### Adapter trimming
    parser.add_argument('-ap','--adaptors',default=False,help="[Optional] A FASTA file containing all adaptors need to be trimmed.")
    ### Parameters for trimmomatic
    parser.add_argument('-tp','--trim_pars',default="",\
        help="[Optional] Arguments and parameters for trimmomatic.")

    ## Annotation
    parser.add_argument('-a','--annotation', default= "", help="[Annotation] Annoteated file for input genome, .GTF format.")
    parser.add_argument('-gn','--genome_name',default="mm10", help="[Annotation] Name of the input genome, default is 'mm10'")
    parser.add_argument('-pl','--plprogram', required=False, help="[Annotation] where is the annotatePeaks.pl")


    args = parser.parse_args()

    main()



    # A filter that only have paritial alignment to chr or AAV 
