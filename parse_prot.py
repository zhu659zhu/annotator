import sys

in_file = open(sys.argv[1], "r")
out_file = open("match_" + sys.argv[1].split("/")[-1], "w")
l = in_file.readline()
while l:
    l1 = in_file.readline()
    # ret_prot = max(l1,l2,l3,key=len).strip()
    ret_prot = l1.strip()
    if len(ret_prot) >= 50:
        out_file.write(l.strip() + "_" + str(len(ret_prot) * 3) + "\n")
        out_file.write(ret_prot + "\n")
    l = in_file.readline()
out_file.close()
in_file.close()
