from subprocess import Popen, PIPE
def run_annotation(bedfile,annfile):

    # bedlines = []

    # with open(bedfile, 'r') as f:
    #     for line in f:
    #         bedlines.append(line.split('\t'))
    # f.close()

    # annlines = {}
    # with open(annfile, 'r') as f:
    #     for line in f:
    #         line_split = line.split('\t')
    #         if line_split[0] not in annlines:
    #             annlines[line_split[0]] = [line_split]   
    #         else:
    #             annlines[line_split[0]].append(line_split)
    # f.close()

    # for i in bedlines:
    #     try:
    #         for j in annlines[i[0]]:
    #             if int(i[1]) > int(j[3]) and int(i[2]) < int(j[4]):
    #                 i.append(j[2])
    #                 i.append(j[1])
    #                 i.append(j[8])
    #     except:
    #         continue

    # f = open(bedfile[:-4] + "_annotated.txt",'w')
    # for line in bedlines:
    #     f.write("\t".join(line))
    # f.close()

    # HOMER VERSION

    filename = bedfile[:-4] + "_annotated.txt"

    arguments = ['perl /Users/mingmingcao/Desktop/Programs/Homer/bin/annotatePeaks.pl {} {} -gtf {} > {}'  
                .format(bedfile,
                        "mm10",
                        annfile,
                        filename)]
    
    process = Popen(args = arguments,
                    shell=True,
                    stdout=PIPE, stderr=PIPE)
    
    stdout, stderr = process.communicate()
    del(stderr)   
    
