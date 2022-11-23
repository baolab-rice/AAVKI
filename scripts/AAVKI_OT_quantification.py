import sys

def run_quantification(filename):

    lines_sam = []
    with open(filename, 'r') as f:
        for line in f: 
            newline = line.strip().split('\t')
            if "chr" in newline[2]:
                lines_sam.append(newline)
    f.close()
    sorted_lines = sorted(lines_sam, key=lambda x: (x[2],int(x[3])))

    bin_size = 1000

    OTs = []
    initial = 0
    i = -1 
    for line in sorted_lines:
        if initial == 0:
            OTs.append([line[2],line[3],1,[line[0]]])
            initial = int(line[3])
            current_chr = line[2]
            i += 1
        else:
            if line[2] == current_chr and initial - bin_size < int(line[3]) < initial + bin_size:
                if line[0] not in OTs[i][3][0]:
                    OTs[i][2] += 1
                    OTs[i][3].append(line[0])
            else: 
                OTs.append([line[2],line[3],1,[line[0]]])
                initial = int(line[3])
                current_chr = line[2]
                i += 1         
    
    f = open(filename[:-4]+".txt",'w')
    f.write("Chr\tPosition\tCount\tRead_ID\n")
    total = 0
    reads = []
    for OT in OTs:
        total += OT[2]
        if OT[2] >5:
            f.write("{}\t{}\t{}\t{}\n".format(OT[0],OT[1],OT[2],OT[3]))
        for read in OT[3]:
            reads.append(read)