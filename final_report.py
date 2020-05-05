
# coding: utf-8

# In[1]:



import os
import time

TE_threshold = 0
protein_evalue_threshold = 1e-50
base_evalue_threshold = 1e-50
protein_length_threshold = 50
base_length_threshold = 150
percent_threshold = 0.05    # translate percentage.
min_translate_length = 60
max_transcript_length = 100000

transcript_dict = {}     # key: transcript_id (e.g. G1.1) value: [gene_id, transcript_id, chr_id, start_pos, end_pos, TE_coverage, cds_similarity, prot_similarity, translate_length, transript_evidence, gene_evidence]
# prot_similarity, cds_similarity: [cs_transcript_id, length, evalue, reverse_sign]

# load bed file and TE_coverage.
os.chdir("./")
in_file = open("KN_9204.bed", "r")
te_in_file = open("result.txt", "r")
l = in_file.readline()
te_l = te_in_file.readline()
while l:
    data = l.strip().split("\t")
    gene_id, transcript_id = data[3].split(";")
    transcript_dict[transcript_id] = [gene_id, transcript_id, data[0], int(data[1]), int(data[2]), float(te_l.strip()), [], [], -1, "", ""]
    l = in_file.readline()
    te_l = te_in_file.readline()
in_file.close()
print "bed file and TE_coverage loaded..."

# load gene evidence
gene_evidence = {}
in_file = open("KN_9204_gene_report.txt", "r")
l = in_file.readline()
l = in_file.readline()
while l:
    data = l.strip().split("\t")
    gene_id = data[0]
    gene_evidence[gene_id] = data[3]
    l = in_file.readline()
in_file.close()
print "gene evidence loaded..."

# load transcript evidence
in_file = open("KN_9204_trans_report.txt", "r")
l = in_file.readline()
l = in_file.readline()
while l:
    data = l.strip().split("\t")
    transcript_id = data[0]
    transcript_dict[transcript_id][9] = data[2]
    transcript_dict[transcript_id][10] = gene_evidence[transcript_dict[transcript_id][0]]
    l = in_file.readline()
in_file.close()
print "transcript evidence loaded..."

# load cds_similarity
in_file = open("cds_origin_KN9204.fasta.outfmt6", "r")
l = in_file.readline()
while l:
    data = l.strip().split("\t")
    transcript_id = data[0].split("_")[0]
#     translate_length = int(data[0].split("_")[1])
    cds_similarity = transcript_dict[transcript_id][6]
#     if len(cds_similarity) == 0:
#         transcript_dict[transcript_id][8] = translate_length
    if len(cds_similarity) == 0 or float(data[10]) < cds_similarity[2]:
        reverse_sign = "+" if int(data[8]) < int(data[9]) else "-"
        cds_similarity = [data[1], int(data[3]), float(data[10]), reverse_sign]
        transcript_dict[transcript_id][6] = cds_similarity
    l = in_file.readline()
in_file.close()
print "cds_similarity loaded..."

# load prot_similarity
in_file = open("match_prot_origin_KN9204.fasta.out", "r")
l = in_file.readline()
while l:
    data = l.strip().split("\t")
    transcript_id = data[0].split("_")[0].replace(">", "")
    prot_similarity = transcript_dict[transcript_id][7]
    cds_similarity = transcript_dict[transcript_id][6]
    if len(prot_similarity) == 0 or float(data[10]) < prot_similarity[2]:
        reverse_sign = "+" if int(data[8]) < int(data[9]) else "-"
        prot_similarity = [data[1], int(data[3]), float(data[10]), reverse_sign]
        transcript_dict[transcript_id][7] = prot_similarity
    l = in_file.readline()
in_file.close()
print "prot_similarity loaded..."

# load transcript length
in_file = open("prot_origin_KN9204.fasta", "r")
l = in_file.readline()
transcript_prot_dict = {}
while l:
    if l[0] == ">":
        transcript_id = l.split(" ")[0].replace(">", "")
        seq = in_file.readline().strip()
        transcript_dict[transcript_id][8] = len(seq) * 3
        transcript_prot_dict[transcript_id] = seq
    l = in_file.readline()
in_file.close()
print "transcript length loaded..."

print transcript_dict["G1.1"]
print transcript_dict["G4.1"]


# In[ ]:


TE_transcript_set = set()

match_gene_chr_set = set()

transcript_set = set()
gene_set = set()
perfect_gene_set = set()
for k, v in transcript_dict.items():
    cds_similarity = v[6]
    prot_similarity = v[7]
    gene_ev = v[10].split(",")
    if (cds_similarity and prot_similarity and cds_similarity[0][7] == v[2][7] and cds_similarity[3] != "-" and cds_similarity[0] == prot_similarity[0] and cds_similarity[2] <= base_evalue_threshold and prot_similarity[2] <= protein_evalue_threshold):
        # 101608
        if "LC" not in cds_similarity[0]:
            match_gene_chr_set.add(cds_similarity[0].split(".")[0])
        perfect_gene_set.add(v[0])
        transcript_set.add(v[1])
        gene_set.add(v[0])
    elif (cds_similarity and prot_similarity and cds_similarity[3] != "-" and cds_similarity[2] <= base_evalue_threshold ** 3 and prot_similarity[2] <= base_evalue_threshold ** 3 and v[5] <= TE_threshold):
        transcript_set.add(v[1])
        gene_set.add(v[0])
    elif (v[9] == "iso_seq" and len(gene_ev) >= 3 and v[5] <= TE_threshold):
        transcript_set.add(v[1])
        gene_set.add(v[0])
    else:
        if v[9] == "iso_seq" and len(gene_ev) >= 2 and v[5] > TE_threshold:
            TE_transcript_set.add(v[1])
        # LC genes
        pass
print len(transcript_set)
print len(gene_set)
print len(perfect_gene_set)

# test on 2019.02.20     need to add v[2] to element of perfect_gene_set !!!!!!!!!!!!!!!!!
match_gene_chr_dict = {}
d4_list = []
for i in match_gene_chr_set:
    chr_id = i[7:9]
    match_gene_chr_dict[chr_id] = match_gene_chr_dict.get(chr_id, 0) + 1
    if chr_id == "4D":
        d4_list.append(i)
print sorted(d4_list)[-50:]
for k in sorted(match_gene_chr_dict.keys()):
    print k, ":", match_gene_chr_dict[k]
exit()




# HC genes: 299594    109081


# for k, v in transcript_dict.items():
#     cds_similarity = v[6]
#     prot_similarity = v[7]
#     gene_ev = v[10].split(",")
#     if (cds_similarity and cds_similarity[2] <= 1e-150 and float(abs(cds_similarity[1] - v[8]))/v[8] < 0.3 and v[0] not in gene_set):
#         # 5376
#         transcript_set.add(v[1])
#         gene_set.add(v[0])
# print len(transcript_set)
# print len(gene_set)



# In[ ]:


# for transcript_id in transcript_set:
#     if transcript_dict[transcript_id][4] - transcript_dict[transcript_id][3] > max_transcript_length:
#         print transcript_id


out_file = open("KN9204_HC_transcript.txt", "w")
gene_unique_isoform_dict = {}
for transcript_id in transcript_set:
    gene_id = transcript_id.split(".")[0]
    isoform_list, unique_isoform_dict = gene_unique_isoform_dict.get(gene_id, [[], {}])
    prot_dict = {}
    # k = transcript_dict[transcript_id][6][0] if transcript_dict[transcript_id][6] else transcript_id  # 被匹配的基因编号
    k = transcript_dict[transcript_id][6][0] + ("_iso" if v[9] == "iso_seq" else "") if transcript_dict[transcript_id][6] else transcript_id  # 被匹配的基因编号
    #k = transcript_dict[transcript_id][6][0] + "_" + str(transcript_dict[transcript_id][6][1]) if transcript_dict[transcript_id][6] else transcript_id  # 被匹配的基因编号
    tmp = unique_isoform_dict.get(k, [])  # old
    v = transcript_dict[transcript_id]  # new
    if v[4] - v[3] > max_transcript_length:
        continue
    if v[8] < min_translate_length:
        continue
    if prot_dict.get(transcript_prot_dict[transcript_id], 0) == 1:
        continue
    
    prot_dict[transcript_prot_dict[transcript_id]] = 1
    unique_isoform_dict[k] = v
#     if not tmp:
#         unique_isoform_dict[k] = v
#     else:
#         if float(v[6][1]) / (v[4] - v[3]) > float(tmp[6][1]) / (tmp[4] - tmp[3]):
#         #if v[6][1] > tmp[6][1] or (v[6][1] == tmp[6][1] and (v[4] - v[3]) < (tmp[4] - tmp[3])):
#             unique_isoform_dict[k] = v
    isoform_list.append(transcript_dict[transcript_id])
    gene_unique_isoform_dict[gene_id] = [isoform_list, unique_isoform_dict]
    
print "loaded..."

for gene_key, gene_value in gene_unique_isoform_dict.items():
    isoform_list, unique_isoform_dict = gene_value
#     if unique_isoform_dict.keys() and ("TraesCS" + unique_isoform_dict.values()[0][2][7] in "".join(unique_isoform_dict.keys())) or ("TraesCSU" in "".join(unique_isoform_dict.keys())):
#         for k, v in unique_isoform_dict.items():
#             if ("TraesCS" + unique_isoform_dict.values()[0][2][7] not in k) and ("TraesCSU" not in k) and (gene_id not in k):
#                 unique_isoform_dict.pop(k)
#     for l in unique_isoform_dict.values():
#         out_file.write(l[1] + "\n")
    for l in sorted(unique_isoform_dict.values(), key=lambda x:int(x[1].split(".")[1])):
        out_file.write(l[1] + "\n")
out_file.close()
print "write finished..."


# In[ ]:


in_file = open("KN9204_HC_transcript.txt", "r")
hc_g_set = set()
hc_t_set = set()
l = in_file.readline()
while l:
    hc_g_set.add(l.split(".")[0])
    hc_t_set.add(l.strip())
    l = in_file.readline()
print len(hc_t_set)
print len(hc_g_set)

tmp_g_set = set()
for transcript_id in hc_t_set:
    prot_similarity = transcript_dict[transcript_id][7]
    cds_similarity = transcript_dict[transcript_id][6]
    cds_evalue = cds_similarity[2] if cds_similarity else 1
    prot_evalue = prot_similarity[2] if prot_similarity else 1
    if cds_evalue > base_evalue_threshold and prot_evalue > protein_evalue_threshold:
        tmp_g_set.add(transcript_id.split(".")[0])
print len(tmp_g_set)


# In[ ]:


# export bed file

id_format = "TraesKN{0}01{1}G{2}0"

out_file_index = open("KN9204_index.txt", "w")
out_file_hc = open("KN9204_HC.bed", "w")
out_file_lc = open("KN9204_LC.bed", "w")
in_file = open("KN_9204.bed", "r")
out_file_TE = open("KN9204_TE.bed", "w")


last_gene_id = ""
index_dict_hc = {}
index_dict_lc = {}
new_gene_index = 0

transcript_set = hc_t_set
print len(transcript_set)
gene_set = set()
for transcript_id in transcript_set:
    gene_id = transcript_id.split(".")[0]
    gene_set.add(gene_id)
print len(gene_set)

for l in in_file:
    data = l.strip().split("\t")
    old_gene_id, old_transcript_id = data[3].split(";")
    
    # export TE gene
    if old_transcript_id in TE_transcript_set:
        out_file_TE.write(l)
    
    if old_gene_id not in gene_set:
        if transcript_dict[old_transcript_id][8] < min_translate_length:
            continue
        confidence_sign = "L"
        if old_gene_id  != last_gene_id:
            last_gene_id = old_gene_id
            index_dict_lc[data[0]] = index_dict_lc.get(data[0], 0) + 1
            new_gene_id = id_format.format(data[0][-2:], confidence_sign, str(index_dict_lc[data[0]]).zfill(4))
            new_gene_index = 1
        else:
            new_gene_index += 1
        data[3] = new_gene_id + ";" + new_gene_id + "." + str(new_gene_index)
        out_file_index.write(old_transcript_id + "\t" + data[3] + "\n")
        out_file_lc.write("\t".join(data) + "\n")
    else:
        if old_transcript_id in transcript_set:
            confidence_sign = "H"
            if old_gene_id  != last_gene_id:
                last_gene_id = old_gene_id
                index_dict_hc[data[0]] = index_dict_hc.get(data[0], 0) + 1
                new_gene_id = id_format.format(data[0][-2:], confidence_sign, str(index_dict_hc[data[0]]).zfill(4))
                new_gene_index = 1
            else:
                new_gene_index += 1
            # data[3] = data[3].replace(old_gene_id, new_gene_id)
            data[3] = new_gene_id + ";" + new_gene_id + "." + str(new_gene_index)
            out_file_index.write(old_transcript_id + "\t" + data[3] + "\n")
            out_file_hc.write("\t".join(data) + "\n")
in_file.close()
out_file_hc.close()
out_file_lc.close()
out_file_TE.close()
out_file_index.close()
print "output..."


# In[ ]:



oldid2newid = {}
newid2oldid = {}
in_file = open("KN9204_index.txt", "r")
l = in_file.readline()
while l:
    data = l.strip().split("\t")
    oldid2newid[data[0]] = data[1]
    newid2oldid[data[1]] = data[0]
    l = in_file.readline()
print "ID index loaded..."


# correct isoform which are independent.
for in_file_path in ["KN9204_HC.bed", "KN9204_LC.bed"]:
    lines = []
    in_file = open(in_file_path, "r")
    for l in in_file:
        lines.append(l.strip().split("\t"))
    lines = sorted(lines, key=lambda x:(x[0], int(x[1])))
    index_dict = {}
    tmp_isoform_index = 1
    new_gene_id = ""
    out_file = open("new_" + in_file_path, "w")
    for i in range(len(lines)):
        line = lines[i]
        if i-1 >= 0 and int(lines[i-1][1]) <= int(line[1]) <= int(lines[i-1][2]):
            tmp_isoform_index += 1
        else:
            index_dict[line[0]] = index_dict.get(line[0], 0) + 1
            tmp_isoform_index = 1
        new_gene_id = id_format.format(line[0][-2:], line[3][11], str(index_dict[line[0]]).zfill(4))
        oldid2newid[newid2oldid[line[3]]] = new_gene_id + ";" + new_gene_id + "." + str(tmp_isoform_index)
        line[3] = new_gene_id + ";" + new_gene_id + "." + str(tmp_isoform_index)
        out_file.write("\t".join(line))
        out_file.write("\n")

    index_file = open(in_file_path + "_index.txt", "w")
    for k,v in oldid2newid.items():
        index_file.write(k+ "\t"+v+"\n")
    print "reindex finished..."
    for i in sorted(index_dict.keys()):
        print i, ": " , index_dict[i]



# In[ ]:


# print transcript_dict["G11822.5"]


# In[ ]:



reverse_gene_set = set()
indel_gene_set = set()
isoseq_gene_set = set()

used_gene_set = set()
used_gene_set.update(perfect_gene_set)

transcript_set = hc_t_set

out_file_reverse = open("spec_gene_reverse.txt", "w")
out_file_indel = open("spec_gene_indel.txt", "w")
out_file_isoseq = open("spec_gene_isoseq.txt", "w")

for transcript_id in transcript_set:
    gene_id = transcript_id.split(".")[0]
    if gene_id in used_gene_set:
        continue
    prot_similarity = transcript_dict[transcript_id][7]
    cds_similarity = transcript_dict[transcript_id][6]
    cds_evalue = cds_similarity[2] if cds_similarity else 1
    prot_evalue = prot_similarity[2] if prot_similarity else 1
    g_id, t_id = oldid2newid[transcript_id].split(";")
    prot_match = prot_similarity[0] if prot_similarity else ""
    cds_match = cds_similarity[0] if cds_similarity else ""
    if cds_similarity and transcript_dict[transcript_id][2][7] == cds_similarity[0][7] and cds_similarity[3] == "-":
        # reverse
        # 711
        # print prot_similarity
        out_file_reverse.write("\t".join([t_id, cds_similarity[0], str(cds_similarity[1]), str(cds_similarity[2])]) + '\n')
        reverse_gene_set.add(gene_id)
        used_gene_set.add(gene_id)
    elif cds_evalue <= base_evalue_threshold and transcript_dict[transcript_id][2][7] == cds_similarity[0][7] and prot_match != cds_match:
        # indel
        # 1363
        if gene_id not in indel_gene_set:
            out_file_indel.write("\t".join([transcript_id, t_id, cds_similarity[0], str(cds_similarity[1]), str(cds_similarity[2])]) + '\n')
            indel_gene_set.add(gene_id)
    elif cds_evalue > base_evalue_threshold and prot_evalue > protein_evalue_threshold:
        # iso-seq
        # 2577
        if "iso_seq" in transcript_dict[transcript_id][10]:
            if gene_id not in isoseq_gene_set:
                out_file_isoseq.write("\t".join([t_id]) + '\n')
                isoseq_gene_set.add(gene_id)
        else:
            print "error"
print len(reverse_gene_set)
print len(indel_gene_set)
print len(isoseq_gene_set)


# In[ ]:


n = 0
for transcript_id in transcript_set:
    if transcript_dict[transcript_id][8] < 60:
        n += 1
print n

# support analyse
gene_support = {}

for transcript_id in transcript_set:
    old_gene_id = transcript_id.split(".")[0]
    new_trans_id = oldid2newid[transcript_id]
    new_gene_id = new_trans_id.split(".")[0]
    gene_support[new_gene_id] = gene_evidence[old_gene_id]

for c in set(gene_support.values()):
    print c, ":", gene_support.values().count(c)



# In[ ]:


# reverse_list = ["1A", "2A", "6A", "1B", "2B", "3B", "5B", "6B", "7B", "1D", "2D", "4D"]
# def process_line(l):
#     chr_length = {
#         "1A": 601718415,
#         "1B": 751739901,
#         "1D": 502780916,
#         "2A": 793025007,
#         "2B": 799391415,
#         "2D": 653434666,
#         "3A": 761766818,
#         "3B": 844935349,
#         "3D": 615516545,
#         "4A": 752191717,
#         "4B": 674887845,
#         "4D": 518212891,
#         "5A": 717244711,
#         "5B": 714813122,
#         "5D": 571824034,
#         "6A": 617560997,
#         "6B": 711369630,
#         "6D": 489148700,
#         "7A": 743621735,
#         "7B": 761607020,
#         "7D": 646941958,
#         "Un": 221792557
#     }
    
#     chr_id = l[0][-2:]
#     if chr_id in reverse_list:
#         l[1], l[2] = chr_length[chr_id] - int(l[2]), chr_length[chr_id] - int(l[1])
#         l[5] = "-" if l[5] == "+" else "+"
#         l[6], l[7] = l[1], l[2]
#         exon_start = l[11].split(",")
#         exon_start = map(int, exon_start)
#         exon_length = l[10].split(",")
#         exon_length = map(int, exon_length)
#         exon_total_length = exon_start[-1] + exon_length[-1]
#         new_exon_start = []
#         for (i1, i2) in zip(exon_start,exon_length):
#             new_exon_start.append(exon_total_length - i1 - i2)
#         new_exon_start.reverse()
#         exon_length.reverse()
#         exon_length = map(str, exon_length)
#         new_exon_start = map(str, new_exon_start)
#         l[10] = ",".join(exon_length)
#         l[11] = ",".join(new_exon_start)
#     else:
#         l[1], l[2] = int(l[1]), int(l[2])
#         l[6], l[7] = l[1], l[2]
        
#     return l
        

# bed_line_list = []

# in_file = open("KN9204_HC.bed", "r")
# for l in in_file:
#     data = l.strip().split("\t")
#     bed_line_list.append(data)

# new_bed_line_list = map(process_line, bed_line_list)
# sorted_new_bed_line_list = []

# for chr_id in ["1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D", "Un"]:
#     tmp_list = filter(lambda x:x[0] == "KN9204." + chr_id, new_bed_line_list)
#     if chr_id in reverse_list:
#         tmp_list.reverse()
#     sorted_new_bed_line_list.extend(tmp_list)

# # sorted_new_bed_line_list = sorted(new_bed_line_list, key=lambda x:(x[0],x[1]))

# id_format = "TraesKN{0}01{1}G{2}0"
# last_gene_id = ""
# index_dict = {}
# new_gene_index = 0

# out_file_index = open("KN9204_new_index.txt", "w")
# out_file = open("new_KN9204_HC.bed", "w")
# for l in sorted_new_bed_line_list:
#     l = map(str, l)
#     old_gene_id = l[3].split(";")[0]
#     if old_gene_id  != last_gene_id:
#         last_gene_id = old_gene_id
#         confidence_sign = old_gene_id[11]
#         index_dict[l[0]] = index_dict.get(l[0], 0) + 1
#         if index_dict.get(l[0], 0) == 1:
#             print old_gene_id
#         new_gene_id = id_format.format(l[0][-2:], confidence_sign, str(index_dict[l[0]]).zfill(4))
#         new_gene_index = 1
#         out_file_index.write(old_gene_id + "\t" + new_gene_id + "\n")
#     else:
#         new_gene_index += 1
#     l[3] = new_gene_id + ";" + new_gene_id + "." + str(new_gene_index)
    
#     out_file.write("\t".join(l))
#     out_file.write("\n")
# out_file.close()
# out_file_index.close()
# print index_dict.items()
# print "output finished..."

