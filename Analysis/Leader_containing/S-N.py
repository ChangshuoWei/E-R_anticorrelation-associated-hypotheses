import os
os.system("samtools view -h Aligned.out.bam > Aligned.out.sam")
fr = open("Aligned.out.sam")
fw = open("junction.sam","w")

for lines in fr :
    if lines.startswith("@"):
        continue
    else :
        if lines.strip().split("\t")[5].count("N") > 0:
            fw.write(lines)
fr.close()
fw.close()
fr = open("junction_filter.sam")
fw = open("junction_read_statistic.txt","w")

S_count = 0
orf3_count = 0
E_count = 0
M_count = 0
orf6_count = 0
orf7a_count = 0
orf8_count = 0
N_count = 0
orf10_count = 0
for lines in fr:
    pos = float(lines.strip().split("\t")[3])
    flag = lines.strip().split("\t")[5]
    if flag.count("S") >0 :
        continue
    if flag.count("D") >0 :
        continue
        temp1 = float(flag.split("M")[0])
        temp2 = float(flag.split("M")[1].split("D")[0])
        temp3 = float(flag.split("M")[1].split("D")[1].split("M")[0])
        site5 = temp1 + temp2 + temp3
        junction = float(flag.split("M")[1].split("D")[1].split("M")[1].split("N")[0])
        site3 = float(pos+site5+junction)
    if flag.count("I") > 0:
        continue
    else:
        site5 = float(flag.split("M")[0])
        junction = float(flag.split("M")[1].split("N")[0])
        site3 = float(pos+site5+junction)
    if pos+site5 < 55 or pos+site5 > 85:
        continue
    if site3 > 21563 - 500 and site3 < 21563 :
        S_count += 1 
    elif site3 > 21563 and  site3 < 25393:
        orf3_count += 1 
    elif  site3 > 25393 and site3 < 26245:
        E_count += 1
    elif site3 > 26245 and  site3 < 26523 :
        M_count += 1 
    elif site3 > 26523 and  site3 < 27202:
        orf6_count += 1 
    elif  site3 > 27202 and  site3 < 27394:
        orf7a_count += 1 
    elif site3 > 27394 and  site3 < 27894:
        orf8_count += 1 
    elif  site3 > 27894 and  site3 < 28274:
        N_count += 1 
    elif site3 > 28274 and site3 < 29558:
        orf10_count +=1

fr.close()
fw.write("S")
fw.write("\t")
fw.write(str(S_count))
fw.write("\n")
fw.write("orf3")
fw.write("\t")
fw.write(str(orf3_count))
fw.write("\n")
fw.write("E")
fw.write("\t")
fw.write(str(E_count))
fw.write("\n")
fw.write("M")
fw.write("\t")
fw.write(str(M_count))
fw.write("\n")
fw.write("orf6")
fw.write("\t")
fw.write(str(orf6_count))
fw.write("\n")
fw.write("orf7a")
fw.write("\t")
fw.write(str(orf7a_count))
fw.write("\n")
fw.write("orf8")
fw.write("\t")
fw.write(str(orf8_count))
fw.write("\n")
fw.write("N")
fw.write("\t")
fw.write(str(N_count))
fw.write("\n")
fw.write("orf10")
fw.write("\t")
fw.write(str(orf10_count))
fw.write("\n")
