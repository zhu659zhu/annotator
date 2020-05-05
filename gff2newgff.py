# coding: utf-8
import sys

in_file = open(sys.argv[1], "r")
#out_file = open("new_" + sys.argv[1].split("/")[-1], "w")
out_file = open("new_KN_9204.gff", "w")
l = in_file.readline()
gene_parent = ""
trans_parent = ""
exon_index = 0
while l:
    if l.strip() == "###":
        pass
    elif l[:2] == "##":
        out_file.write(l)
    else:
        data = l.strip().split("\t")
        name = data[8].split(";")
        if data[2] == "gene":
            if name[2].endswith(".1"):
                out_file.write("###\n")
                gene_parent = name[2].split(".")[0]
                trans_parent = name[2]
                data[8] = "ID=" + gene_parent + ";Name=" + gene_parent
                out_file.write("\t".join(data))
                out_file.write("\n")
                data[8] = "ID=" + trans_parent + ";Parent=" + gene_parent + ";Name=" + trans_parent
                data[2] = "mRNA"
                out_file.write("\t".join(data))
                out_file.write("\n")
            else:
                exon_index = 0
                trans_parent = name[2]
                data[8] = "ID=" + trans_parent + ";Parent=" + gene_parent + ";Name=" + trans_parent
                data[2] = "mRNA"
                out_file.write("\t".join(data))
                out_file.write("\n")
        elif data[2] == "exon":
            exon_index += 1
            data[8] = "ID=" + trans_parent + ".exon" + str(exon_index) + ";Parent=" + trans_parent
            out_file.write("\t".join(data))
            out_file.write("\n")
    l = in_file.readline()

out_file.close()
in_file.close()
