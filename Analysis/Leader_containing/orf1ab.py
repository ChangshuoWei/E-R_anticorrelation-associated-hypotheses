import os 
fr = open("leader.container.sam")
fw = open("orf1ab.sam","w")

for lines in fr:
    if lines.split("\t")[5].count("N"):
        continue
    else:
        fw.write(lines)
fr.close()
fw.close()
os.system("wc -l orf1ab.sam")
