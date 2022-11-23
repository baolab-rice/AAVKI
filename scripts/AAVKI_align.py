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
