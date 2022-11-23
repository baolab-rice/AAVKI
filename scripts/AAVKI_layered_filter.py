import re

def filter_sam():
    pass

def read_sam(sam):

    reads_dict = {}
    with open(sam, 'r') as f:
        for line in f: 
            if '@' in line.split('\t')[0]:
                continue
            else:
                newline = line.split('\t')
                read_ID = newline[0]
                if read_ID not in reads_dict:
                    reads_dict[read_ID] = [line]
                else:
                    reads_dict[read_ID].append(line)
    f.close()

    return reads_dict

def filter_AAV_CRISPR(dictname):
    
    AAV_CRISPR_dict = {}
    nAAV_CRISPR_dict = {}

    for title,reads in dictname.items():
        refs = []
        for read in reads:
            refs.append(read.split('\t')[2])
        if "AAV_CRISPR" in refs:
            AAV_CRISPR_dict[title] = reads
        else:
            nAAV_CRISPR_dict[title] = reads
    
    return AAV_CRISPR_dict,nAAV_CRISPR_dict

def filter_AAV_donor(sam,HR):

    threshold = HR + 10

    reads_AAVR = {}
    reads_genome = {}

    for title,reads in sam.items():
        for read in reads: 
            cigars = []
            if read.split('\t')[2] == "AAV_Donor":
                cigars.append(read.split('\t')[5])
            else:
                if "SA:" in read:
                    SA = read.split('\t')[15]
                    if "AAV_Donor" in SA:
                        cigars.append(re.split(":|,|;",SA)[re.split(":|,|;",SA).index("AAV_Donor") + 3])
                if "XA:" in read:
                    try:
                        XA = read.split('\t')[16]
                    except:
                        try:
                            XA = read.split('\t')[15]
                        except:
                            continue
                    if "AAV_Donor" in XA:
                        cigars.append(re.split(":|,|;",XA)[re.split(":|,|;",XA).index("AAV_Donor") + 2])
            if cigars != []:
                # print(cigars)
                for cigar in cigars:
                    p = re.compile("[A-Z]")
                    start = 0
                    for match in p.finditer(cigar):
                        length = cigar[start:match.start()]
                        start = match.start() + 1
                        if match.group() == "M":
                            if int(length) > threshold and title not in reads_AAVR:
                                reads_AAVR[title] = reads

    for title,reads in sam.items():
        if title not in reads_AAVR:
            reads_genome[title] = reads

    return reads_AAVR,reads_genome

def group_ON_and_OT(genome_dict,site,HR):

    ON_chr = site.split(':')[0]
    ON_site = int(site.split(':')[1])
     
    reads_chrON = {}
    reads_chrOT = {}

    for title,reads in genome_dict.items():
        for read in reads:
            newline = read.split('\t')
            Chr = newline[2]
            mstart = newline[3]
            if Chr == ON_chr and ON_site - 50 < int(mstart) < ON_site + 50:
                cigar = newline[5]
                p = re.compile("[A-Z]")
                start = 0
                for match in p.finditer(cigar):  
                    length = cigar[start:match.start()]
                    start = match.start() + 1
                    if match.group() == "M":
                        if int(length) > HR + 10: 
                            reads_chrON[title] = reads         
                        else:
                            reads_chrOT[title] = reads
            else:
                reads_chrOT[title] = genome_dict[title]

    return reads_chrON,reads_chrOT

def group_by_ITR_containing(dictname):

    reads_ITR = {}
    reads_nITR = {}

    for title,reads in dictname.items():
        refs = []
        for read in reads:
            refs.append(read.split('\t')[2])
        if "AAV_ITR" in refs:
            reads_ITR[title] = reads
        else:
            reads_nITR[title] = reads
    
    return reads_ITR,reads_nITR

def write_output(filename,dictname):

    f = open(filename,'w')
    for ID,reads in dictname.items():
        for read in reads:
            f.write(read)
    f.close()

def run_layered_filters(sam,HR,site):

    reads_dict = read_sam(sam)
    AAV_CRISPR_dict,nAAV_CRISPR_dict = filter_AAV_CRISPR(reads_dict)
    print("1")
    AAV_donor_dict,genome_dict = filter_AAV_donor(nAAV_CRISPR_dict,HR)
    print("2")
    ON_dict,OT_dict = group_ON_and_OT(genome_dict,site,HR)
    print("3")

    ON_ITRs,ON_nITR = group_by_ITR_containing(ON_dict)
    OT_ITRs,OT_nITR = group_by_ITR_containing(OT_dict)
    print("4")

    # Outputs
    fb = sam[:-4]

    ## recombination w/ AAV_CRISPR
    AAV_CRISPR_file = fb + "_AAV_CRISPR.sam"
    write_output(AAV_CRISPR_file,AAV_CRISPR_dict)

    ## recombination w/ AAV_donor 
    AAV_donor_file = fb + "_AAV_donor_nCRISPR.sam"
    write_output(AAV_donor_file,AAV_donor_dict)

    ## ON & OT
    ON_HDR_file = fb + "_ON_HDR.sam"
    ON_NHEJ_file = fb + "_ON_NHEJ.sam" 
    OT_HDR_file = fb + "_OT_HDR.sam" 
    OT_NHEJ_file = fb + "_OT_NHEJ.sam" 

    write_output(ON_HDR_file,ON_nITR)
    write_output(ON_NHEJ_file,ON_ITRs)
    write_output(OT_HDR_file,OT_nITR)
    write_output(OT_NHEJ_file,OT_ITRs)

    return OT_HDR_file,OT_NHEJ_file