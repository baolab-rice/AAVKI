from subprocess import Popen, PIPE

def run_align(reads,genome):

    filename = reads[:-6] + ".sam"

    arguments = ['bwa mem {} {} > {}'  
                .format(genome,
                        reads,
                        filename)]
    
    process = Popen(args = arguments,
                    shell=True,
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    del(stderr)   

    return filename

def samtobed(samfile):

    bamfile = samfile[:-4] + ".bam"

    arguments = ['samtools view -bS {} > {}'  
                .format(samfile,
                        bamfile)]
    
    process = Popen(args = arguments,
                    shell=True,
                    stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    del(stderr)

    bedfile = bamfile[:-4] + ".bed"

    arguments = ["bedtools bamtobed -i {} > {}"
                 .format(bamfile,
                         bedfile)]
    
    process = Popen(args = arguments,
                    shell=True,
                    stdout=PIPE, stderr=PIPE)
    
    stdout, stderr = process.communicate()
    del(stderr)  

    return bedfile



