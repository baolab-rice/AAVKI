import sys
import re

def is_polynucleotide(seq, threshold):

    bases = set("ACGT")
    n_bases = len(seq)
    n_polynuc = max([seq.count(b) for b in bases])
    return (n_polynuc / n_bases) >= threshold


def run_quantification(filename,site):

    ON_chr = site.split(':')[0]
    ON_site = int(site.split(':')[1])

    lines_sam = []

    lines_sam_dict = {}
    AAVreITR_list = []
    filtered_reads = []
    reads_with_AAV_CRISPR = []

    with open(filename, 'r') as f:
        for line in f: 
            newline = line.strip().split('\t')
            if newline[2] == "AAV_CRISPR":
                if newline[0] not in reads_with_AAV_CRISPR:
                    reads_with_AAV_CRISPR.append(newline[0])
                    continue
            ID = newline[0]
            if "AAV_Donor" in newline[2]:
                continue    
            if "AAV_ITR" in newline[2]:
                AAVreITR_list.append(line)
                continue
            elif ON_chr in newline[2]:
                if ON_site - 1200 < int(newline[3]) and int(newline[3]) < ON_site + 1200:
                    continue
                # v4 previous score of 60 
                elif int(newline[4]) >= 0:
                    lines_sam.append(newline)
                    if ID not in lines_sam_dict:
                        lines_sam_dict[ID] = [newline]
                    else:
                        lines_sam_dict[ID].append(newline)
            # The "*" filter also has in the layered filter
            elif newline[2] != "*" and int(newline[4]) >= 0:
                lines_sam.append(newline)
                if ID not in lines_sam_dict:
                    lines_sam_dict[ID] = [newline]
                else:
                    lines_sam_dict[ID].append(newline)
            else:
                filtered_reads.append(line)

    # print(len(filtered_reads))

    # for ID, value in lines_sam_dict.items():
    #     print(ID,value)

    f.close()
    sorted_lines = sorted(lines_sam, key=lambda x: (x[2],int(x[3])))

    # for line in sorted_lines:
    #     print(line)

    bin_size = 1000
    threshold = 0.60

    OTs = []
    initial = 0
    i = -1 
    p = re.compile("[A-Z]")

    quantified_reads = []
    multi_reads = []
    low_length_match = []
    
    temp_check = []


    for line in sorted_lines:

        ID = line[0]
        if ID not in temp_check:
            temp_check.append(ID)
        seq_full = line[9]
        cigar = line[5]
        start = 0
        start_cigar = 0

        # didn't check I or D in 2 Ms
        for match in p.finditer(cigar):
            if match.group() == "S": 
                start += int(cigar[start:match.start()])    
                start_cigar = match.start() + 1
            elif match.group() == "M":
                end = start + int(cigar[start_cigar:match.start()]) - 1
                length = int(cigar[start_cigar:match.start()])
                seq = seq_full[start:end]
                # print(seq_full,start,cigar,end,seq)
                break
            else:
                start_cigar = match.start() + 1

        ## Filter of matching length
        if length < 50:
            # print(ID, length)
            if ID not in low_length_match:
                low_length_match.append(ID)
            continue
        else:
            if ID not in quantified_reads:
                quantified_reads.append(ID)
            else:
                multi_reads.append(ID)

            if is_polynucleotide(seq,threshold) == True:
                is_poly = "Poly"
            else:
                is_poly = "Non-poly"
            
            if initial == 0:
                OTs.append([line[2],line[3],1,[line[0]],[line[5]],[seq],is_poly,[]])
                initial = int(line[3])
                current_chr = line[2]
                i += 1
            else:
                if line[2] == current_chr:
                    if initial - bin_size < int(line[3]) and int(line[3]) < initial + bin_size:
                        if line[0] not in OTs[i][3][0]:
                            OTs[i][2] += 1
                            OTs[i][3].append(line[0])
                            OTs[i][4].append(line[5])
                            OTs[i][5].append(seq)
                    else:
                        OTs.append([line[2],line[3],1,[line[0]],[line[5]],[seq],is_poly,[]])
                        initial = int(line[3])
                        current_chr = line[2]
                        i += 1     
                else: 
                    OTs.append([line[2],line[3],1,[line[0]],[line[5]],[seq],is_poly,[]])
                    initial = int(line[3])
                    current_chr = line[2]
                    i += 1   

    # some reads may contain multiplt hits can counted twice
    for ID in quantified_reads:
        if ID in low_length_match:
            low_length_match.remove(ID)

    # Output the read with short length matching to geonme
    f = open(filename[:-4]+"_filtered_short_match.txt",'w')
    for ID, reads in lines_sam_dict.items():
        if ID in low_length_match:
            f.write("\n".join("\t".join(map(str,l)) for l in reads))
            f.write("\n")
    f.close()

    # label reads with AAV_CRISPR
    for ID in reads_with_AAV_CRISPR:
        for OT in OTs:
            if ID in OT[3]:
                OT[7].append((ID,"w/AAV_CRISPR"))

    f = open(filename[:-4]+".txt",'w')
    f.write("Chr\tPosition\tCount\tIf_poly\tRead_ID\tCigars\tMatched_seq\twAAV_CRISPR\n")
    total = 0
    reads = []
    k = 1 
    for OT in OTs:
        total += OT[2]
        if OT[2] > 0:
            if OT[6] == "Non-poly":
                title = "BINNED_SEQ_" + str(k)
                f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(OT[0],OT[1],OT[2],OT[6],OT[3],OT[4],OT[5],OT[7],title))
                k += 1
            else:
                f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(OT[0],OT[1],OT[2],OT[6],OT[3],OT[4],OT[5],OT[7]))
        for read in OT[3]:
            reads.append(read)
    print("Reads number: {}".format(len(quantified_reads)))
    print("OT Events: {}".format(total))
    print("Reads filtered of short matching length < 50 bp: {}".format(len(low_length_match)))
    print("Multi_hits_reads: {}".format(multi_reads))
    f.close()

    f = open(filename[:-4]+".fasta",'w')
    k = 1
    for OT in OTs:
        if OT[6] == "Non-poly":
            seq = sorted(OT[5], reverse=True)[0]
            title = "BINNED_SEQ_" + str(k)
            f.write(">{}\n".format(title))
            f.write(seq)
            f.write("\n")
            k += 1
    f.close()

    return filename[:-4]+".fasta"


def run_ON_quantification(filename,site):
    ON_chr = site.split(':')[0]
    ON_site = int(site.split(':')[1])

    lines_sam = []
    lines_sam_dict = {}
    reads_with_AAV_CRISPR = []
    with open(filename, 'r') as f:
        for line in f: 
            newline = line.strip().split('\t')
            ID = newline[0]
            # AAV_CRISPR checking
            if newline[2] == "AAV_CRISPR":
                if newline[0] not in reads_with_AAV_CRISPR:
                    reads_with_AAV_CRISPR.append(newline[0])
                    
            if  ON_chr in newline[2]:                        
                if ON_site - 1200 > int(newline[3]) or int(newline[3]) > ON_site + 1200:
                    continue
                else:
                    lines_sam.append(newline)
                if ID not in lines_sam_dict:
                    lines_sam_dict[ID] = [newline]
                else:
                    lines_sam_dict[ID].append(newline)
    f.close()

    sorted_lines = sorted(lines_sam, key=lambda x: (x[2],int(x[3])))

    bin_size = 10   
    OTs = []
    initial = 0
    i = -1 
    p = re.compile("[A-Z]")

    quantified_reads = []
    

    for line in sorted_lines:

        if line[2] == "AAV_CRISPR":
            print("yes")
            if line[0] not in reads_with_AAV_CRISPR:
                reads_with_AAV_CRISPR.append(line[0])
                continue
        else:
            if line[0] not in quantified_reads:
                quantified_reads.append(line[0])
            else:
                continue

            if initial == 0:
                OTs.append([line[2],line[3],1,[line[0]],[line[5]],[]])
                initial = int(line[3])
                current_chr = line[2]
                i += 1
            else:
                if line[2] == current_chr:
                    if initial - bin_size < int(line[3]) and int(line[3]) < initial + bin_size:
                        if line[0] not in OTs[i][3][0]:
                            OTs[i][2] += 1
                            OTs[i][3].append(line[0])
                            OTs[i][4].append(line[5])
                    else:
                        OTs.append([line[2],line[3],1,[line[0]],[line[5]],[]])
                        initial = int(line[3])
                        current_chr = line[2]
                        i += 1     
                else: 
                    OTs.append([line[2],line[3],1,[line[0]],[line[5]],[]])
                    initial = int(line[3])
                    current_chr = line[2]
                    i += 1   

    # label reads with AAV_CRISPR
    for ID in reads_with_AAV_CRISPR:
        for OT in OTs:
            if ID in OT[3]:
                OT[5].append((ID,"w/AAV_CRISPR"))

    f = open(filename[:-4]+".txt",'w')
    f.write("Chr\tPosition\tCount\tRead_ID\tCigars\twAAV_CRISPR\n")
    total = 0
    reads = []
    for OT in OTs:
        total += OT[2]
        if OT[2] > 0:
            f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(OT[0],OT[1],OT[2],OT[3],OT[4],OT[5]))
        for read in OT[3]:
            reads.append(read)
    print("Reads number: {}".format(len(lines_sam_dict)))
    print("ON Events: {}".format(total))
    f.close()

    # Feature 1 need to consider: Tracking each raed ID
    # Feature 2 need to consider: mapping size filter (e.g > 30bp)
    # Feature 3 ployA mapping trimming
