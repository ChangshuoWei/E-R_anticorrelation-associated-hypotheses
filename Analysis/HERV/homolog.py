fr = open("query_result_Hsap")
fw = open("homolog_id_Hsap.txt",'w')

judge_seq = ''
for lines in fr:
    query_seq = lines.strip().split("\t")[0]
    if query_seq != judge_seq:
        judge_seq = query_seq
    else:
        continue
    refer_seq = lines.strip().split("\t")[1]
    coverage = lines.strip().split("\t")[2]
    if float(coverage) < 85:
        continue
    query_chr = query_seq.split(".")[1]
    refer_chr = refer_seq.split(".")[1]
    if refer_chr == "chr2a" or refer_chr == "chr2b":
        if query_chr== "chr2" :
            fw.write(lines)
    elif query_chr == refer_chr:
        fw.write(lines)
fr.close()
fw.close()
import os
fr = open("homolog_id_Hsap.txt")
fw = open("homolog_info_Hsap.txt",'w')
fw.write("query_id")
fw.write("\t")
fw.write("refer_id")
fw.write("\t")
fw.write("Identities")
fw.write("\t")
fw.write("cover_length")
fw.write("\t")
fw.write("different_aa")
fw.write("\t")
fw.write("gap")
fw.write("\t")
fw.write("query_start")
fw.write("\t")
fw.write("query_end")
fw.write("\t")
fw.write("refer_start")
fw.write("\t")
fw.write("refer_end")
fw.write("\t")
fw.write("E-value")
fw.write("\t")
fw.write("score")
fw.write("\t")
fw.write("query_seq")
fw.write("\t")
fw.write("refer_seq")
fw.write("\t")
fw.write("query_length")
fw.write("\t")
fw.write("refer_length")
fw.write("\n")
for lines in fr:
    query_id = lines.strip().split("\t")[0]
    refer_id = lines.strip().split("\t")[1]
    query_temp =  os.popen("grep '"+query_id+"' -A 1 Hsap38.fa")
    refer_temp =  os.popen("grep '"+refer_id+"' -A 1 Ptro214.fa")
    query_seq = query_temp.read().split("\n")[1]
    refer_seq = refer_temp.read().split("\n")[1]
    fw.write(lines.strip())
    fw.write("\t")
    fw.write(query_seq)
    fw.write("\t")
    fw.write(refer_seq)
    fw.write("\t")
    fw.write(str(len(query_seq)))
    fw.write("\t")
    fw.write(str(len(refer_seq)))
    fw.write("\n")
    
fr.close()
fw.close()

fr = open("query_result_Ptro")
fw = open("homolog_id_Ptro.txt",'w')

judge_seq = ''
for lines in fr:
    query_seq = lines.strip().split("\t")[0]
    if query_seq != judge_seq:
        judge_seq = query_seq
    else:
        continue
    refer_seq = lines.strip().split("\t")[1]
    coverage = lines.strip().split("\t")[2]
    if float(coverage) < 85:
        continue
    query_chr = query_seq.split(".")[1]
    refer_chr = refer_seq.split(".")[1]
    if refer_chr == "chr2":
        if query_chr == "chr2a" or query_chr == "chr2b":
            fw.write(lines)
    elif query_chr == refer_chr:
        fw.write(lines)
fr.close()
fw.close()
import os
fr = open("homolog_id_Ptro.txt")
fw = open("homolog_info_Ptro.txt",'w')
fw.write("query_id")
fw.write("\t")
fw.write("refer_id")
fw.write("\t")
fw.write("Identities")
fw.write("\t")
fw.write("cover_length")
fw.write("\t")
fw.write("different_aa")
fw.write("\t")
fw.write("gap")
fw.write("\t")
fw.write("query_start")
fw.write("\t")
fw.write("query_end")
fw.write("\t")
fw.write("refer_start")
fw.write("\t")
fw.write("refer_end")
fw.write("\t")
fw.write("E-value")
fw.write("\t")
fw.write("score")
fw.write("\t")
fw.write("query_seq")
fw.write("\t")
fw.write("refer_seq")
fw.write("\t")
fw.write("query_length")
fw.write("\t")
fw.write("refer_length")
fw.write("\n")
for lines in fr:
    query_id = lines.strip().split("\t")[0]
    refer_id = lines.strip().split("\t")[1]
    query_temp =  os.popen("grep '"+query_id+"' -A 1 Hsap38.fa")
    refer_temp =  os.popen("grep '"+refer_id+"' -A 1 Ptro214.fa")
    query_seq = query_temp.read().split("\n")[1]
    refer_seq = refer_temp.read().split("\n")[1]
    fw.write(lines.strip())
    fw.write("\t")
    fw.write(query_seq)
    fw.write("\t")
    fw.write(refer_seq)
    fw.write("\t")
    fw.write(str(len(query_seq)))
    fw.write("\t")
    fw.write(str(len(refer_seq)))
    fw.write("\n")
    
fr.close()
fw.close()
