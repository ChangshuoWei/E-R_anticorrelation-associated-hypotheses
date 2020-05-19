import os
os.system("makeblastdb -in refer_seq.fa -dbtype nucl -out refer_seq")
os.system("blastn -query gisaid_0512_7310.fa -db refer_seq -outfmt 3 -out query_results")
fr = open('query_results')
fw = open('check_var.txt','w')
count = 0
sample_name = "sample"
for lines in fr:
    if lines.startswith('Query='):
        sample_name = lines.strip().split("=")[1].strip()
        fw.write(sample_name)
        fw.write('\n')
        count = 0
    elif lines.startswith("Query_"):
        query_seq = lines.split()[2]
    elif lines.startswith('0'):
        if lines[1] != ' ':
            continue
        refer_seq = lines.split()[2]
        refer_start = lines.split()[1]
        refer_end = lines.split()[3]
        for i in range(len(refer_seq)):
            if refer_seq[i] == '.':
                continue
            else :
                count += 1
                refer_position= int(float(refer_start)+float(i))
                if query_seq[i].upper() in ['A','T','C','G'] and refer_seq[i] != '-':
                    fw.write(str(refer_position))
                    fw.write('\t')
                    fw.write(refer_seq[i])
                    fw.write('\t')
                    fw.write(query_seq[i])
                    fw.write('\n')
    else :
        continue
fr.close()
fw.close()
fr = open("check_var.txt")
fw = open("abstract.txt","w")

for lines in fr:
    if lines.startswith("hCoV"):
        continue
    else:
        fw.write(lines.upper())
fr.close()
fw.close()
