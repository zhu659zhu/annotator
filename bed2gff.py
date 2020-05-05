# coding: utf-8
import sys

gene_dict = {}

#for j in range(1,2):
#    for i in ["A"]:
for j in range(1, 8):
    for i in ["A", "B", "D"]:
         with open("/home/sczhuhd/data_wheat_kn9204/KN9204_ref/final_KN9204/KN9204." + str(j)+i + ".fasta","r") as in_file:
             in_file.readline()
             l = in_file.readline()
             temp = ""
             while l :
                 temp += l.strip()
                 l = in_file.readline()
             gene_dict[str(j)+i] = temp
             print "finished loading Chr" + str(j) + i
             print len(temp)

j = "U"
i = "n"
with open("/home/sczhuhd/data_wheat_kn9204/KN9204_ref/final_KN9204/KN9204." + str(j)+i + ".fasta","r") as in_file:
     in_file.readline()
     l = in_file.readline()
     temp = ""
     while l :
         temp += l.strip()
         l = in_file.readline()
     gene_dict[str(j)+i] = temp
     print "finished loading Chr" + str(j) + i
     print len(temp)
print gene_dict.keys()

# tips: gff file start from 1, including the end. bed file start from 0, excluding the end.
def get_sequence(scaffold_id, start, end):
    # for gff file.
    start = int(start)
    end = int(end)
    temp = gene_dict[scaffold_id]
    rsp = temp[start:end]
    return rsp

def get_protein(seq):
    seq = seq.upper()
    codonTable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'.', 'TAG':'.',
    'TGC':'C', 'TGT':'C', 'TGA':'.', 'TGG':'W',
    }
    ret = []
    ret_seq = []
    ret_loc = []
    for i in range(0, 3):
        tmp_prot = ""
        tmp_seq = ""
        tmp_loc = 0
        start = 0
        for loc in range(i, len(seq) - 2, 3):
            codon = seq[loc:loc+3]
            if codon == "ATG" and start == 0:
                tmp_loc = loc
                start = 1
            if start == 1:
                tmp_seq += codon
                tmp_prot = tmp_prot + codonTable.get(codon, "*")
            if start == 1 and codonTable.get(codon, "") == ".":
                ret.append(tmp_prot)
                ret_seq.append(tmp_seq)
                ret_loc.append(tmp_loc)
                start = 0
                tmp_prot = ""
                tmp_seq = ""
        ret.append(tmp_prot)
        ret_seq.append(tmp_seq)
        ret_loc.append(tmp_loc)
    cds = max(ret_seq, key=len)
    index = ret_seq.index(cds)
    prot = ret[index]
    loc = ret_loc[index]
    return cds, prot, loc

def get_reverse(seq):
    def get_nuc(x):
        x = x.upper()
        if x == "A":
            return "T"
        elif x == "T":
            return "A"
        elif x == "C":
            return "G"
        elif x == "G":
            return "C"
        elif x == "N":
            return "N"
    seq = seq.upper()
    ret = ""
    for i in seq:
        ret = get_nuc(i) + ret
    return ret


def get_cds_range(list_a,list_b,start,length,sep='+'):
    seq_range = []
    first_loc = list_a[0]
    for i in range(list_a[-1]+list_b[-1]-first_loc):
        seq_range.append(0)
    for (exon_start, exon_length) in zip(list_a, list_b):
        for i in range(exon_start, exon_start + exon_length):
            seq_range[i - first_loc] = 1
    if sep == "+":
        start_index = 0
        length_index = length
        for i in range(len(seq_range)):
            if start_index >= start and seq_range[i] == 1 and length_index > 0:
                length_index -= 1
                seq_range[i - first_loc] = 2
            if start_index < start and seq_range[i] == 1:
                start_index += 1
            
    if sep == "-":
        start_index = 0
        length_index = length
        for i in range(len(seq_range)-1, -1, -1):
            if start_index >= start and seq_range[i] == 1 and length_index > 0:
                length_index -= 1
                seq_range[i - first_loc] = 2
            if start_index < start and seq_range[i] == 1:
                start_index += 1
    # print "".join(map(str, [x for x in seq_range]))
    
    ans_l = []
    ans = []
    ans_r = []
    range_sign = 0
    tmp_s = 0
    now_base = ""
    for index in range(len(seq_range)):
        tmp_sign = seq_range[index]
        if now_base != tmp_sign:
            if now_base != 0:
                if range_sign == 0:
                    range_sign = 1
                elif range_sign == 1:
                    ans_l.append([tmp_s, index])
                elif range_sign == 2:
                    ans.append([tmp_s, index])
                elif range_sign == 3:
                    ans_r.append([tmp_s, index])
            if tmp_sign == 2:
                range_sign = 2
            elif ans:
                range_sign = 3
            elif tmp_sign == 1:
                range_sign = 1
            tmp_s = index
            now_base = tmp_sign
    if now_base == 2:
        ans.append([tmp_s, index + 1])
    else:
        ans_r.append([tmp_s, index + 1])
    if sep == "-":
        ans_l, ans_r = ans_r, ans_l
        
    xr_ans_l, xr_ans, xr_ans_r = get_cds_range_xr(list_a,list_b,start,length,sep)
    if xr_ans_l!=ans_l or xr_ans != ans or xr_ans_r != ans_r:
        print list_a,list_b,start,length,sep
        print ans_l, ans, ans_r
        print xr_ans_l, xr_ans, xr_ans_r
    
    return ans_l, ans, ans_r

    
def get_cds_range_xr(list_a,list_b,start,length,sep='+'):
    def dup_filter(x):
        if x[0] != x[1]:
            return True
        else:
            return False

    intervals = []
    for i in range(len(list_a)):
        intervals.append([list_a[i], list_a[i]+list_b[i]])

    begin = 0
    ans_start = []
    ans = []
    ans_end = []


    if sep == '-':
        intervals = intervals[::-1]
    i = 0
    ok = 1
    yes = 1
    while i < len(intervals):
        begin += intervals[i][1] - intervals[i][0]
        if begin < start:
            ans_start.append(intervals[i])
        elif begin >= start and begin < start + length:
            if ok == 1:
                pre_sum = sum([d[1]-d[0] for d in ans_start])
                if sep == '+':
                    ans_start.append((intervals[i][0], intervals[i][0]+start-pre_sum))
                    ans.append((intervals[i][0]+start-pre_sum, intervals[i][1]))
                else:
                    ans_start.append((intervals[i][1]-(start-pre_sum),intervals[i][1]))
                    ans.append((intervals[i][0], intervals[i][1]-(start-pre_sum)))
                ok = 0
            else:
                ans.append(intervals[i])

        else:
            if ok == 1 and yes == 1: 
                if sep == '+':
                    pre_sum1 = sum([d[1]-d[0] for d in ans_start])
                    ans_start.append((intervals[i][0], intervals[i][0]+start-pre_sum1))
                    ok = 0
                    # print(pre_sum1)
                    ans.append((start-pre_sum1+intervals[i][0], start-pre_sum1+length+intervals[i][0]))
                    yes = 0
                    pre_sum2 = sum([d[1]-d[0] for d in ans_start + ans])
                    # print(pre_sum2)
                    ans_end.append((length+start-pre_sum1+intervals[i][0], intervals[i][1]))
                else:
                    pre_sum1 = sum([d[1]-d[0] for d in ans_start])
                    ans_start.append((intervals[i][1]-(start-pre_sum1),intervals[i][1]))
                    ans.append((intervals[i][1]-(start+length-pre_sum1),intervals[i][1]-(start-pre_sum1)))
                    # print(pre_sum1)
                    ok = 0
                    pre_sum2 = sum([d[1]-d[0] for d in ans_start + ans])
                    # print(pre_sum2)
                    yes = 0
                    ans_end.append((intervals[i][0], intervals[i][1]-(start-pre_sum1+length)))
                
            elif yes == 1:
                if sep == '+':
                    pre_sum = sum([d[1]-d[0] for d in ans_start + ans])
                    ans.append((intervals[i][0], intervals[i][0]+start+length-pre_sum))
                    ans_end.append((intervals[i][0]-pre_sum+start+length, intervals[i][1]))
                else:
                    pre_sum = sum([d[1]-d[0] for d in ans_start + ans])
                    ans.append((intervals[i][1]-(start+length-pre_sum),intervals[i][1]))
                    ans_end.append((intervals[i][0], intervals[i][1]-(start+length-pre_sum)))
                    
                yes = 0
            else:
                ans_end.append(intervals[i])
        i += 1
        
    for i in range(len(ans_start)):
        ans_start[i] = list(ans_start[i])
    for i in range(len(ans)):
        ans[i] = list(ans[i]) 
    for i in range(len(ans_end)):
        ans_end[i] = list(ans_end[i])        
    ans_start = list(filter(dup_filter, ans_start))
    ans = list(filter(dup_filter, ans))
    ans_end = list(filter(dup_filter, ans_end))
    ans_start.sort(key=lambda x:x[0])
    ans.sort(key=lambda x:x[0])
    ans_end.sort(key=lambda x:x[0])
    return ans_start,ans, ans_end


gff_header_dict = {}  # max and min location.

out_file_name = "IGDB_v1.0_" + sys.argv[1].split("/")[-1].replace(".bed", "") +"_201811"

in_file = open(sys.argv[1], "r")
out_file = open(out_file_name + ".gff3", "w")
out_file_cds = open(out_file_name + "_cds" + ".fasta", "w")
out_file_prot = open(out_file_name+ "_pep" + ".fasta", "w")
out_file_transcript = open(out_file_name+ "_transcripts" + ".fasta", "w")
out_file_gtf = open(out_file_name + ".gtf", "w")
l = in_file.readline()
gene_parent = ""
trans_parent = ""
exon_index = 0
#out_file.write("##gff-version 3\n")

version_str = "IGDB_v1.0_201811"

while l:
    data = l.strip().split("\t")
    gene_id, transcript_id = data[3].split(";")
    chr_id = data[0][-2:]
    gff_header_dict[chr_id+"_min"] = min(gff_header_dict.get(chr_id+"_min", 999999999), int(data[1]))
    gff_header_dict[chr_id+"_max"] = max(gff_header_dict.get(chr_id+"_max", 0), int(data[2]))
    if transcript_id.endswith(".1"):
        if gene_parent != "":
            out_file.write("###\n")
        gene_parent = gene_id
        out_file.write("{}\t{}\t{}\t{}\t{}\t.\t{}\t.\t{}\n".format(data[0],version_str , "gene", str(int(data[1]) + 1), str(int(data[2])), data[5], "ID="+gene_id+";Name="+gene_id))
    trans_parent = transcript_id
    out_file_transcript.write(">" + transcript_id + " gene=" + gene_parent + "\n")
    out_file.write("{}\t{}\t{}\t{}\t{}\t.\t{}\t.\t{}\n".format(data[0], version_str, "mRNA", str(int(data[1]) + 1), str(int(data[2])), data[5], "ID="+transcript_id+";Parent=" + gene_id+";Name="+transcript_id))
    seq = ""
    exon_length = data[10].split(",")
    exon_length = map(int ,exon_length)
    exon_start = data[11].split(",")
    exon_start = map(int, exon_start)
    chr_id = data[0][-2:]
    for i in range(int(data[9])):
        seq += get_sequence(chr_id, exon_start[i]+int(data[1]), exon_start[i]+int(data[1])+exon_length[i])
        out_file.write("{}\t{}\t{}\t{}\t{}\t.\t{}\t.\t{}\n".format(data[0], version_str, "exon", str(int(data[1]) + 1 + exon_start[i]), str(exon_start[i]+int(data[1])+exon_length[i]), data[5], "ID="+transcript_id+".exon"+str(i+1)+";Parent=" + transcript_id))
        out_file_gtf.write("{}\t{}\t{}\t{}\t{}\t.\t{}\t.\t{}\n".format(data[0], version_str, "exon", str(int(data[1]) + 1 + exon_start[i]), str(exon_start[i]+int(data[1])+exon_length[i]), data[5], "transcript_id \"" + transcript_id + "\"; gene_id \"" + gene_parent + "\"; gene_name \"" + gene_parent + "\""))
    if data[5] == "-":
        seq = get_reverse(seq)
    out_file_transcript.write(seq + "\n")
    cds_seq, prot_seq, start_loc = get_protein(seq)
    five_prime_utr, tmp_cds_list, three_prime_utr = get_cds_range(exon_start, exon_length, start_loc, len(cds_seq), data[5])
    
    seg_list = []
    for cds in tmp_cds_list:
        seg_list.append(str(cds[0]) + "-" + str(cds[1]))
    seg_str = ",".join(seg_list)
    #out_file_cds.write(">"+ trans_parent + " gene=" + gene_parent + " loc:chr" + chr_id + "\n")
    out_file_cds.write(">%s gene=%s loc:chr%s(%s)%s-%s segs:%s\n" % (trans_parent, gene_parent, chr_id, data[5], str(int(data[1]) + 1), str(int(data[2])), seg_str))
    out_file_cds.write(cds_seq + "\n")
    out_file_prot.write(">"+ trans_parent + " gene=" + gene_parent + "\n")
    out_file_prot.write(prot_seq + "\n")
    for i, cds in enumerate(tmp_cds_list):
        out_file_gtf.write("{}\t{}\t{}\t{}\t{}\t.\t{}\t.\t{}\n".format(data[0], version_str, "CDS", str(int(data[1]) + 1 + cds[0]), str(cds[1] +int(data[1])), data[5], "transcript_id \""+transcript_id + "\"; gene_id \"" + gene_parent + "\"; gene_name \"" + gene_parent + "\""))
    if data[5] == "+":
        for i, utr in enumerate(five_prime_utr):
            out_file.write("{}\t{}\t{}\t{}\t{}\t.\t{}\t.\t{}\n".format(data[0], version_str, "five_prime_UTR", str(int(data[1]) + 1 + utr[0]), str(utr[1] + int(data[1])), data[5], "ID="+transcript_id+".utr5p"+str(i+1)+";Parent=" + transcript_id))
        for i, cds in enumerate(tmp_cds_list):
            out_file.write("{}\t{}\t{}\t{}\t{}\t.\t{}\t.\t{}\n".format(data[0], version_str, "CDS", str(int(data[1]) + 1 + cds[0]), str(cds[1] +int(data[1])), data[5], "ID="+transcript_id+".cds"+str(i+1)+";Parent=" + transcript_id))
        for i, utr in enumerate(three_prime_utr):
            out_file.write("{}\t{}\t{}\t{}\t{}\t.\t{}\t.\t{}\n".format(data[0], version_str, "three_prime_UTR", str(int(data[1]) + 1 + utr[0]), str(utr[1] + int(data[1])), data[5], "ID="+transcript_id+".utr3p"+str(i+1)+";Parent=" + transcript_id))
    else:
        for i, utr in enumerate(three_prime_utr):
            out_file.write("{}\t{}\t{}\t{}\t{}\t.\t{}\t.\t{}\n".format(data[0], version_str, "three_prime_UTR", str(int(data[1]) + 1 + utr[0]), str(utr[1] + int(data[1])), data[5], "ID="+transcript_id+".utr3p"+str(i+1)+";Parent=" + transcript_id))
        for i, cds in enumerate(tmp_cds_list):
            out_file.write("{}\t{}\t{}\t{}\t{}\t.\t{}\t.\t{}\n".format(data[0], version_str, "CDS", str(int(data[1]) + 1 + cds[0]), str(cds[1] +int(data[1])), data[5], "ID="+transcript_id+".cds"+str(i+1)+";Parent=" + transcript_id))
        for i, utr in enumerate(five_prime_utr):
            out_file.write("{}\t{}\t{}\t{}\t{}\t.\t{}\t.\t{}\n".format(data[0], version_str, "five_prime_UTR", str(int(data[1]) + 1 + utr[0]), str(utr[1] + int(data[1])), data[5], "ID="+transcript_id+".utr5p"+str(i+1)+";Parent=" + transcript_id))

    l = in_file.readline()
out_file.close()
in_file.close()

# add gff3 header
out_file = open(out_file_name + ".gff3", "r+")
content = out_file.read()
out_file.seek(0, 0)
out_file.write("##gff-version 3\n")
for k in ["1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D", "Un"]:
    out_file.write("##sequence-region   KN9204.%s %s %s\n" % (k, gff_header_dict[k + "_min"], gff_header_dict[k + "_max"]))
out_file.write(content)
out_file.close()

