from subprocess import Popen, PIPE
from Bio import SeqIO

def trim_reads(read1,read2,adaptors,pars):

    threads = 4

    arguments = ['trimmomatic PE -threads {} {} {} {} {} {} {} ILLUMINACLIP:{}/:2:30:10 {}'  
                .format(
                        threads,
                        read1,
                        read2,
                        read1_trim,
                        read1_untrim,
                        read2_trim,
                        read2_untrim,
                        adaptors,
                        pars)]
    
    process = Popen(args = arguments,
                    shell=True,
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    del(stderr)   

def merge_reads(read1,read2):

    ## FLASH 1.2.11 was used

    arguments = ['flash -M600 {} {} 2>&1 | tee flash.log'  
                .format(read1,
                        read2)]
    
    process = Popen(args = arguments,
                    shell=True,
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()

    del(stderr) 

# def cat_paired_reads(read1,read2):

#     read1 = SeqIO.to_dict(SeqIO.parse(read1,"fastq"))
#     read2 = SeqIO.to_dict(SeqIO.parse(read2,"fastq"))

#     filename = path_to_sample + "/" + file_ID + "_unmerged_filtered.fasta"
#     f = open(filename,'w')

#     for read in read1:
#         if bc in read1[read].seq:
#             f.write(">" + read + "\n")
#             f.write(str(read1[read].seq) + str(read2[read].seq) + "\n")
    
#     return filename

def filter_single_read(read,bc):

    def to_dict_remove_dups(sequences):
        return {record.id: record for record in sequences}

    
    reads = to_dict_remove_dups(SeqIO.parse(read, "fastq"))

    filename = path_to_sample + file_ID + "_filtered.fasta"

    f = open(filename,'w')

    for read in reads:
        if bc in reads[read].seq:
            f.write(">" + read + "\n")
            f.write(str(reads[read].seq) + "\n")  

    # If no filter 
    # for read in reads:
    #     f.write(">" + read + "\n")
    #     f.write(str(reads[read].seq) + "\n")  
    
    return filename

# def filter_non_specific_PCR_product(bc):

#     read1 = "./out.notCombined_1.fastq"
#     read2 = "./out.notCombined_2.fastq"
#     read  = "./out.extendedFrags.fastq"

#     unmerged = cat_paired_reads(read1,read2)
#     merged = filter_single_read(read)

#     arguments = ['cat {} >> {}'  
#                 .format(unmerged,
#                         merged)]
    
#     process = Popen(args = arguments,
#                     shell=True,
#                     stdout=PIPE, stderr=PIPE)
#     stdout, stderr = process.communicate()
#     del(stderr) 

#     return merged

def run_QC_and_merge(read2,adaptors,pars,bc):
    # read1_trim,read1_untrim,

    global read2_trim,read2_untrim,path_to_sample,file_ID

    # if read1[-3:] == ".gz":
    #     read1_trim = read1[:-9] + "_trimed.fastq"
    #     read1_untrim = read1[:-9] + "_untrimed.fastq"
    # else:
    #     read1_trim = read1[:-6] + "_trimed.fastq"
    #     read1_untrim = read1[:-6] + "_untrimed.fastq"        

    if read2[-3:] == ".gz":
        read2_trim = read2[:-9] + "_trimed.fastq"
        read2_untrim = read2[:-9] + "_untrimed.fastq"
    else:
        read2_trim = read2[:-6] + "_trimed.fastq"
        read2_untrim = read2[:-6] + "_untrimed.fastq"    

    path_to_sample = read2[:(len(read2) - len(read2.split("/")[-1]))]

    file_ID = read2.split("/")[-1]
    file_ID = file_ID.split(".fastq")[0]

    # trim_reads(read1,read2,adaptors,pars)

    # merged_reads = filter_single_read(read2_trim,bc)
    merged_reads = filter_single_read(read2,bc)

    # merge_reads(read1_trim,read2_trim)
    # merged_reads = filter_non_specific_PCR_product(bc)

    return merged_reads
