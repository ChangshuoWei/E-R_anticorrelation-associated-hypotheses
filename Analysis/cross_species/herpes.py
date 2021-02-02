import os
import pandas as pd
import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
fr = open("name.txt")
name_list = []
for lines in fr:
    name_list.append(lines.strip())
fr.close()
fw = open("beta_cds.fasta",'w')
for name in name_list:
    fr_temp = open(name+"/"+name+"_cds.fasta")
    fr_temp.readline()
    seq = fr_temp.readline()
    fw.write(">"+name)
    fw.write("\n")
    fw.write(seq)
    fr_temp.close()
fw.close()
fr = open("beta.txt")
fw = open("beta_best.txt",'w')
flag = " "
for lines in fr:
    ref_name = lines.strip().split("\t")[0]
    if ref_name == flag:
        continue
    flag = ref_name
    fw.write(lines)
fr.close()
fw.close()
fr = open("Panine.txt")
fw = open("Panine_best.txt",'w')
flag = " "
for lines in fr:
    ref_name = lines.strip().split("\t")[0]
    if ref_name == flag:
        continue
    flag = ref_name
    fw.write(lines)
fr.close()
fw.close()
fr = open("Panine_best.txt")
reciprical_dict = {}
for lines in fr:
    name = lines.strip().split("\t")[0]
    reciprical_dict[name] = lines.strip().split("\t")[1]
fr.close()
fr = open("beta_best.txt")
fw = open("reciprical_best_hit.txt",'w')
for lines in fr:
    name1 = lines.strip().split("\t")[0]
    name2 = lines.strip().split("\t")[1]
    if name1 == reciprical_dict[name2]:
        fw.write(lines)
fr.close()
fw.close()
fr=open('Panine_cds.fasta', 'r')
fw=open('Panine_cds.fa', 'w')
seq={}
for line in fr:
    if line.startswith('>'):    #判断字符串是否以‘>开始’
        name=line.strip()    #以空格为分隔符，并取序列为0的项。
        seq[name]=''
    else:
        seq[name]+=line.replace('\n', '')
fr.close()                           

for i in seq.keys():
    fw.write(i)
    fw.write('\n')
    fw.write(seq[i])
    fw.write('\n')
fw.close()
fr = open("reciprical_best_hit.txt")
for lines in fr:
    name1 = lines.strip().split("\t")[0]
    name2 = lines.strip().split("\t")[1]
    seq1 = os.popen("grep '"+name1+"' -A 1 beta_cds.fa")
    seq1 = seq1.read()
    seq2 = os.popen("grep '"+name2+"' -A 1 Panine_cds.fa")
    seq2 = seq2.read()
    file_name_pre = name1
    fw = open("cds_seq/"+file_name_pre+"_cds.fasta",'w')
    id1 = seq1.split("\n")[0].strip()
    id2 = seq2.split("\n")[0].strip()
    beta_seq = Seq(seq1.split("\n")[1].strip())
    Panine_seq = Seq(seq2.split("\n")[1].strip())
    beta_seq_cds = str(beta_seq)
    fw.write(id1)
    fw.write("\n")
    fw.write(beta_seq_cds)
    fw.write("\n")
    if id2.count("complement") >0:
        Panine_seq_cds = str(Panine_seq.reverse_complement())
        fw.write(id2)
        fw.write("\n")
        fw.write(Panine_seq_cds)
        fw.write("\n")
    else:
        Panine_seq_cds = str(Panine_seq)
        fw.write(id2)
        fw.write("\n")
        fw.write(Panine_seq_cds)
        fw.write("\n")
    fw.close()
    fw = open("protein_seq/"+file_name_pre+"_pro.fasta",'w')

    beta_seq_pro = str(beta_seq.translate())
    fw.write(id1)
    fw.write("\n")
    fw.write(beta_seq_pro)
    fw.write("\n")
    if id2.count("complement") >0:
        Panine_seq_pro = str(Panine_seq.reverse_complement().translate())
        fw.write(id2)
        fw.write("\n")
        fw.write(Panine_seq_pro)
        fw.write("\n")
    else:
        Panine_seq_pro = str(Panine_seq.translate())
        fw.write(id2)
        fw.write("\n")
        fw.write(Panine_seq_pro)
        fw.write("\n")
    fw.close()
fr.close()
fw = open("PanineVSbeta.txt",'w')
fw.write("protein")
fw.write("\t")
fw.write("dnds")
fw.write("\t")
fw.write("dn")
fw.write("\t")
fw.write("ds")
fw.write("\t")
fw.write("ts/tv")
fw.write("\n")
fr = open("gene_id.txt")
name_list = []
for lines in fr:
    name_list.append(lines.strip())
fr.close()
for name_ID in name_list:
    cds_file_path = "cds_seq/"+name_ID+"_cds.fasta"
    pro_file_path = "protein_seq/"+name_ID+"_pro.fasta"
    pro_align_file_path = "protein_seq/"+name_ID+"_pro_align.fasta"
    cds_align_file_path = "cds_seq/"+name_ID+"_cds_align.fasta"
    cds_seq = SeqIO.parse(cds_file_path,"fasta")
    pro_seq = SeqIO.parse(pro_file_path,"fasta")
    i = 0
    for seq_record in cds_seq:
        if i == 0:
            beta_cds = seq_record
        else:
            Panine_cds = seq_record
        i += 1
    i = 0
    for seq_record in pro_seq:
        if i == 0:
            beta_pro = seq_record
        else:
            Panine_pro = seq_record
        i += 1
    pro_cds_dict_beta = {}
    beta_pro_seq = str(beta_pro.seq)
    beta_cds_seq = str(beta_cds.seq)
    for i in range(len(beta_pro_seq)):
        name = beta_pro_seq[i]+str(i)
        pro_cds_dict_beta[name] = beta_cds_seq[3*i:3*(i+1)]
    pro_cds_dict_Panine = {}
    Panine_pro_seq = str(Panine_pro.seq)
    Panine_cds_seq = str(Panine_cds.seq)
    for i in range(len(Panine_pro_seq)):
        name = Panine_pro_seq[i]+str(i)
        pro_cds_dict_Panine[name] = Panine_cds_seq[3*i:3*(i+1)]
    os.system("muscle -in "+pro_file_path+" -out "+pro_align_file_path)
    protein_seq = SeqIO.parse(pro_align_file_path,"fasta")
    i = 0
    for seq_record in protein_seq:
        if i == 0:
            beta_pro_align = seq_record
        else:
            Panine_pro_align = seq_record
        i += 1
    i = 0
    beta_cds_align = ''
    for j in range(len(beta_pro_align.seq)):
        if str(beta_pro_align.seq)[j] != '-':
            name = str(beta_pro_align.seq)[j]+str(i)
            beta_cds_align += pro_cds_dict_beta[name]
            i += 1
        elif str(beta_pro_align.seq)[j] == '-':
            beta_cds_align += '---'
    i = 0
    Panine_cds_align = ''
    for j in range(len(Panine_pro_align.seq)):
        if str(Panine_pro_align.seq)[j] != '-':
            name = str(Panine_pro_align.seq)[j]+str(i)
            Panine_cds_align += pro_cds_dict_Panine[name]
            i += 1
        elif str(Panine_pro_align.seq)[j] == '-':
            Panine_cds_align += '---'
    fw_temp = open("temp.seq",'w')
    fw_temp.write("  2  "+str(len(beta_cds_align)))
    fw_temp.write("\n")
    fw_temp.write("beta")
    fw_temp.write("\n")
    fw_temp.write(beta_cds_align)
    fw_temp.write("\n")
    fw_temp.write("Panine")
    fw_temp.write("\n")
    fw_temp.write(Panine_cds_align)
    fw_temp.close()
    os.system("codeml temp.ctl")
    fr_temp = open("temp.Null.result")
    for lines in fr_temp:
        if lines.startswith("kappa"):
            tstv = lines.strip().split("=")[1].strip()
        elif lines.startswith("omega"):
            dnds = lines.strip().split("=")[1].strip()
        elif lines.startswith("tree length for dN"):
            dn = lines.strip().split(":")[1].strip()
        elif lines.startswith("tree length for dS"):
            ds = lines.strip().split(":")[1].strip()
    fr_temp.close()
    fw.write(name_ID)
    fw.write("\t")
    fw.write(str(dnds))
    fw.write("\t")
    fw.write(str(dn))
    fw.write("\t")
    fw.write(str(ds))
    fw.write("\t")
    fw.write(str(tstv))
    fw.write("\n")
fw.close()
fr = open("name.txt")
fw = open("gene_id.txt",'w')
for lines in fr:
    fw.write(lines.strip().split("_cds")[0])
    fw.write("\n")
fr.close()
fw.close()
