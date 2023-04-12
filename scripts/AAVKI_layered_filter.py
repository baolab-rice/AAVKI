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

def filter_AAV_CRISPR(dictname,site,HR):

    ON_chr = site.split(':')[0]
    ON_site = int(site.split(':')[1])

    AAV_CRISPR_dict = {}
    nAAV_CRISPR_dict = {}

    for title,reads in dictname.items():
        is_AAV_CRISPR = False
        for read in reads:
            read_split = read.split('\t')
            if read_split[2] == "AAV_CRISPR":
                is_AAV_CRISPR = True
                for read_ in reads:
                    read_split_ = read_.split('\t')
                    if read_split_[2] == ON_chr:
                        cigar = read_split_[5]
                        p = re.compile("[A-Z]")
                        start = 0
                        for match in p.finditer(cigar):  
                            length = cigar[start:match.start()]
                            start = match.start() + 1
                            if match.group() == "M":
                                if int(length) > HR + 10: 
                                    nAAV_CRISPR_dict[title] = reads 
                    elif "chr" in read_split_[2] and read_split_[2] != ON_chr:
                        nAAV_CRISPR_dict[title] = reads
        if is_AAV_CRISPR == False:
            nAAV_CRISPR_dict[title] = reads
    
    for title,reads in dictname.items():
        if title not in nAAV_CRISPR_dict.keys():
            AAV_CRISPR_dict[title] = reads

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
                # SA tag for marking chimeric reads
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

def filter_nonsense_reads(genome_dict,site,HR):

    ON_chr = site.split(':')[0]
    ON_site = int(site.split(':')[1])

    to_delete = []

    filtered_genome_dict = genome_dict.copy()

    for title,reads in genome_dict.items():
        mapping_chr = []
        cigars = []
        for read in reads:
            newline = read.split('\t')
            Chr = newline[2]
            ### may over filter (doing)
            # if Chr != "AAV_Donor" and Chr != ON_chr:
            #     if "AAV_Donor" in read or ON_chr in read:
            #         Chr = "AAV_Donor"
            cigar = newline[5]
            mapping_chr.append(Chr)
            cigars.append(cigar)
        # 1st filter for mapping only to HR
        ## may consider underestimate 
        if set(mapping_chr).difference({"AAV_Donor",ON_chr}) == set():
            p = re.compile("[A-Z]")
            marker = 0
            for cigar in cigars:
                start = 0
                for match in p.finditer(cigar):
                    length = cigar[start:match.start()]
                    start = match.start() + 1
                    if match.group() == "M":
                        # On target
                        if int(length) > HR + 10:
                            marker = 1
            if marker == 0:
                to_delete.append(title)
        # 2nd filter for not mapping to HR
        elif "AAV_Donor" not in mapping_chr and ON_chr not in mapping_chr:
            to_delete.append(title)
    
    for title in to_delete:
        del filtered_genome_dict[title]

    return filtered_genome_dict


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
            # Need to consider the values
            if Chr == ON_chr:
                marker = 0
                if ON_site - 100 < int(mstart) and int(mstart) < ON_site + 100:
                    cigar = newline[5]
                    p = re.compile("[A-Z]")
                    start = 0
                    for match in p.finditer(cigar):  
                        length = cigar[start:match.start()]
                        start = match.start() + 1
                        if match.group() == "M":
                            if int(length) > HR + 10: 
                                reads_chrON[title] = reads 
                                marker = 1 
                # Add another filter to v4
                if marker == 0:
                    # when chr9, could filter with AAV_ITR "AAV_Donor" not in read and "AAV_ITR" not in read and newline[2] != "*" and 
                    if abs(int(mstart) - ON_site) > 1200:
                        reads_chrOT[title] = reads
            elif "chr" in Chr:
                reads_chrOT[title] = genome_dict[title]

    keys = genome_dict.keys() - reads_chrOT.keys() - reads_chrON.keys()

    other_AAVre_dict = {}
    for key in keys:
        other_AAVre_dict[key] = genome_dict[key]

    return reads_chrON,reads_chrOT, other_AAVre_dict

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
    AAV_CRISPR_dict,nAAV_CRISPR_dict = filter_AAV_CRISPR(reads_dict,site,HR)
    AAV_donor_dict,genome_dict = filter_AAV_donor(nAAV_CRISPR_dict,HR)

    # remove reads that non-sense (be careful)
    filter_genome_dict = filter_nonsense_reads(genome_dict,site,HR)

    ON_dict,OT_dict,other_AAVre_dict = group_ON_and_OT(filter_genome_dict,site,HR)

    final_genome_dict = dict(list(ON_dict.items()) + list(OT_dict.items()))

    ON_ITRs,ON_nITR = group_by_ITR_containing(ON_dict)
    OT_ITRs,OT_nITR = group_by_ITR_containing(OT_dict)

    # Outputs
    fb = sam[:-4]

    ## recombination w/ AAV_CRISPR
    AAV_CRISPR_file = fb + "_AAV_CRISPR.sam"
    write_output(AAV_CRISPR_file,AAV_CRISPR_dict)

    ## recombination w/ AAV_donor 
    AAV_donor_file = fb + "_AAV_donor_nCRISPR.sam"
    write_output(AAV_donor_file,AAV_donor_dict)

    ## other AAV recombination
    other_AAVre_file = fb + "_other_AAV_recombination.sam"
    write_output(other_AAVre_file,other_AAVre_dict)

    ## w/ genome 
    AAV_genome_file = fb + "_on_genome.sam"
    AAV_genome_filtered_file = fb + "_on_genome_filtered.sam"

    write_output(AAV_genome_file, genome_dict)
    write_output(AAV_genome_filtered_file, filter_genome_dict)

    ## ON & OT
    ON_HDR_file = fb + "_ON_HDR.sam"
    ON_NHEJ_file = fb + "_ON_NHEJ.sam" 
    OT_HDR_file = fb + "_OT_HDR.sam" 
    OT_NHEJ_file = fb + "_OT_NHEJ.sam" 

    write_output(ON_HDR_file,ON_nITR)
    write_output(ON_NHEJ_file,ON_ITRs)
    write_output(OT_HDR_file,OT_nITR)
    write_output(OT_NHEJ_file,OT_ITRs)

    ## Reads number stats
    print("Reads total: {}".format(len(reads_dict)))
    print("Reads with long matching to AAV_CRISPR: {}".format(len(AAV_CRISPR_dict)))
    print("Reads in nAAV_CRISPR: {}".format(len(nAAV_CRISPR_dict)))
    print("Reads with long matching to AAV_Donor: {}".format(len(AAV_donor_dict)))
    print("Filtered nonsense reads (mapping to nowhere): {}".format(len(genome_dict)-len(filter_genome_dict)))
    print("Reads are considered as AAV recombination: {}".format(len(other_AAVre_dict)))
    print("Reads sent for on genome analysis: {}".format(len(final_genome_dict)))
    print("Reads for OT_HDR quantification: {}".format(len(OT_nITR)))
    print("Reads for OT_NHEJ quantification: {}".format(len(OT_ITRs)))
    print("Reads for ON_HDR quantification: {}".format(len(ON_nITR)))
    print("Reads for ON_NHEJ quantification: {}".format(len(ON_ITRs)))

    return OT_HDR_file,OT_NHEJ_file,ON_HDR_file,ON_NHEJ_file