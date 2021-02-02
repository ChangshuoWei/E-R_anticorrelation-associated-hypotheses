import os
import pandas as pd
import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
fr = open("vaccinia.txt")
fw = open("vaccinia_best.txt",'w')
flag = " "
for lines in fr:
    ref_name = lines.strip().split("\t")[0]
    if ref_name == flag:
        continue
    flag = ref_name
    fw.write(lines)
fr.close()
fw.close()
fr = open("shunkpox.txt")
fw = open("shunkpox_best.txt",'w')
flag = " "
for lines in fr:
    ref_name = lines.strip().split("\t")[0]
    if ref_name == flag:
        continue
    flag = ref_name
    fw.write(lines)
fr.close()
fw.close()
fr = open("shunkpox_best.txt")
reciprical_dict = {}
for lines in fr:
    name = lines.strip().split("\t")[0]
    reciprical_dict[name] = lines.strip().split("\t")[1]
fr.close()
fr = open("vaccinia_best.txt")
fw = open("reciprical_best_hit.txt",'w')
for lines in fr:
    name1 = lines.strip().split("\t")[0]
    name2 = lines.strip().split("\t")[1]
    if name1 == reciprical_dict[name2]:
        fw.write(lines)
fr.close()
fw.close()
fr = open("shunkpox_cds.fasta")
fw = open("shunkpox_cds.fa",'w')
seq = {}
for lines in fr:
    if lines.startswith(">"):
        name = lines.strip()
        seq[name] = ''
    else:
        seq[name] += lines.strip()
fr.close()
for key in seq.keys():
    fw.write(key)
    fw.write("\n")
    fw.write(seq[key])
    fw.write("\n")
fw.close()
fr = open("reciprical_best_hit.txt")
for lines in fr:
    name1 = lines.strip().split("\t")[0]
    name2 = lines.strip().split("\t")[1]
    seq1 = os.popen("grep '"+name1+"' -A 1 vaccinia_cds.fa")
    seq1 = seq1.read()
    seq2 = os.popen("grep '"+name2+"' -A 1 shunkpox_cds.fa")
    seq2 = seq2.read()
    file_name_pre = name1.split(".1_")[1]
    fw = open("cds_seq/"+file_name_pre+"_cds.fasta",'w')
    id1 = seq1.split("\n")[0].strip()
    id2 = seq2.split("\n")[0].strip()
    vaccinia_seq = Seq(seq1.split("\n")[1].strip())
    shunkpox_seq = Seq(seq2.split("\n")[1].strip())
    vaccinia_seq_cds = str(vaccinia_seq)
    fw.write(id1)
    fw.write("\n")
    fw.write(vaccinia_seq_cds)
    fw.write("\n")
    if id2.count("complement") >0:
        shunkpox_seq_cds = str(shunkpox_seq.reverse_complement())
        fw.write(id2)
        fw.write("\n")
        fw.write(shunkpox_seq_cds)
        fw.write("\n")
    else:
        shunkpox_seq_cds = str(shunkpox_seq)
        fw.write(id2)
        fw.write("\n")
        fw.write(shunkpox_seq_cds)
        fw.write("\n")
    fw.close()
    fw = open("protein_seq/"+file_name_pre+"_pro.fasta",'w')

    vaccinia_seq_pro = str(vaccinia_seq.translate())
    fw.write(id1)
    fw.write("\n")
    fw.write(vaccinia_seq_pro)
    fw.write("\n")
    if id2.count("complement") >0:
        shunkpox_seq_pro = str(shunkpox_seq.reverse_complement().translate())
        fw.write(id2)
        fw.write("\n")
        fw.write(shunkpox_seq_pro)
        fw.write("\n")
    else:
        shunkpox_seq_pro = str(shunkpox_seq.translate())
        fw.write(id2)
        fw.write("\n")
        fw.write(shunkpox_seq_pro)
        fw.write("\n")
    fw.close()
fr.close()
fw = open("shunkpoxVSvaccinia.txt",'w')
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
            vaccinia_cds = seq_record
        else:
            shunkpox_cds = seq_record
        i += 1
    i = 0
    for seq_record in pro_seq:
        if i == 0:
            vaccinia_pro = seq_record
        else:
            shunkpox_pro = seq_record
        i += 1
    pro_cds_dict_vaccinia = {}
    vaccinia_pro_seq = str(vaccinia_pro.seq)
    vaccinia_cds_seq = str(vaccinia_cds.seq)
    for i in range(len(vaccinia_pro_seq)):
        name = vaccinia_pro_seq[i]+str(i)
        pro_cds_dict_vaccinia[name] = vaccinia_cds_seq[3*i:3*(i+1)]
    pro_cds_dict_shunkpox = {}
    shunkpox_pro_seq = str(shunkpox_pro.seq)
    shunkpox_cds_seq = str(shunkpox_cds.seq)
    for i in range(len(shunkpox_pro_seq)):
        name = shunkpox_pro_seq[i]+str(i)
        pro_cds_dict_shunkpox[name] = shunkpox_cds_seq[3*i:3*(i+1)]
    os.system("muscle -in "+pro_file_path+" -out "+pro_align_file_path)
    protein_seq = SeqIO.parse(pro_align_file_path,"fasta")
    i = 0
    for seq_record in protein_seq:
        if i == 0:
            vaccinia_pro_align = seq_record
        else:
            shunkpox_pro_align = seq_record
        i += 1
    i = 0
    vaccinia_cds_align = ''
    for j in range(len(vaccinia_pro_align.seq)):
        if str(vaccinia_pro_align.seq)[j] != '-':
            name = str(vaccinia_pro_align.seq)[j]+str(i)
            vaccinia_cds_align += pro_cds_dict_vaccinia[name]
            i += 1
        elif str(vaccinia_pro_align.seq)[j] == '-':
            vaccinia_cds_align += '---'
    i = 0
    shunkpox_cds_align = ''
    for j in range(len(shunkpox_pro_align.seq)):
        if str(shunkpox_pro_align.seq)[j] != '-':
            name = str(shunkpox_pro_align.seq)[j]+str(i)
            shunkpox_cds_align += pro_cds_dict_shunkpox[name]
            i += 1
        elif str(shunkpox_pro_align.seq)[j] == '-':
            shunkpox_cds_align += '---'
    fw_temp = open("temp.seq",'w')
    fw_temp.write("  2  "+str(len(vaccinia_cds_align)))
    fw_temp.write("\n")
    fw_temp.write("vaccinia")
    fw_temp.write("\n")
    fw_temp.write(vaccinia_cds_align)
    fw_temp.write("\n")
    fw_temp.write("shunkpox")
    fw_temp.write("\n")
    fw_temp.write(shunkpox_cds_align)
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
fr = open("shunkpoxVSvaccinia.txt")
fw = open("shunkpoxVSvaccinia_process.txt",'w')
for lines in fr:
    fw.write(lines.replace('_',''))
fr.close()
fw.close()
fr = open("gene_id.txt")
name_list = []
for lines in fr:
    name_list.append(lines.strip())
fr.close()
for name in name_list:
    file_path = "protein_seq/"+name+"_pro.fasta"
    fw= open(file_path,'w')
    temp = os.popen("cat ../vaccinia/protein_seq/"+name+"_pro.fasta")
    fw.write(temp.read())
    temp = os.popen("tail -2 ../vaccinia_smallpox/protein_seq/"+name+"_pro.fasta")
    fw.write(temp.read())
    temp = os.popen("tail -2 ../vaccinia_monkey/protein_seq/"+name+"_pro.fasta")
    fw.write(temp.read())
    temp = os.popen("tail -2 ../vaccinia_skunk/protein_seq/"+name+"_pro.fasta")
    fw.write(temp.read())
    fw.close()
fr = open("gene_id.txt")
name_list = []
for lines in fr:
    name_list.append(lines.strip())
fr.close()
for name in name_list:
    file_path = "cds_seq/"+name+"_cds.fasta"
    fw= open(file_path,'w')
    temp = os.popen("cat ../vaccinia/cds_seq/"+name+"_cds.fasta")
    fw.write(temp.read())
    temp = os.popen("tail -2 ../vaccinia_smallpox/cds_seq/"+name+"_cds.fasta")
    fw.write(temp.read())
    temp = os.popen("tail -2 ../vaccinia_monkey/cds_seq/"+name+"_cds.fasta")
    fw.write(temp.read())
    temp = os.popen("tail -2 ../vaccinia_skunk/cds_seq/"+name+"_cds.fasta")
    fw.write(temp.read())
    fw.close()
fw = open("vaccinia_tree.txt",'w')
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
for name_ID in name_list:
    cds_file_path = "cds_seq/"+name_ID+"_cds.fasta"
    pro_file_path = "protein_seq/"+name_ID+"_pro.fasta"
    pro_align_file_path = "protein_seq/"+name_ID+"_pro_align.fasta"
    os.system("muscle -in "+pro_file_path+" -out "+pro_align_file_path)
    protein_seq = SeqIO.parse(pro_align_file_path,"fasta")
    cds_seq = SeqIO.parse(cds_file_path,"fasta")
    pro_seq = SeqIO.parse(pro_file_path,"fasta")
    i = 0
    for seq_record in cds_seq:
        if i == 0:
            vaccinia_cds = seq_record
        elif i == 1:
            cowpox_cds = seq_record
        elif i == 2:
            smallpox_cds = seq_record
        elif i == 3:
            monkeypox_cds = seq_record
        elif i == 4:
            skunkpox_cds = seq_record
        i += 1
    i = 0
    for seq_record in pro_seq:
        if i == 0:
            vaccinia_pro = seq_record
        elif i == 1:
            cowpox_pro = seq_record
        elif i == 2:
            smallpox_pro = seq_record
        elif i == 3:
            monkeypox_pro = seq_record
        elif i == 4:
            skunkpox_pro = seq_record
        i += 1
    pro_cds_dict_vaccinia = {}
    vaccinia_pro_seq = str(vaccinia_pro.seq)
    vaccinia_cds_seq = str(vaccinia_cds.seq)
    for i in range(len(vaccinia_pro_seq)):
        name = vaccinia_pro_seq[i]+str(i)
        pro_cds_dict_vaccinia[name] = vaccinia_cds_seq[3*i:3*(i+1)]
    pro_cds_dict_cowpox = {}
    cowpox_pro_seq = str(cowpox_pro.seq)
    cowpox_cds_seq = str(cowpox_cds.seq)
    for i in range(len(cowpox_pro_seq)):
        name = cowpox_pro_seq[i]+str(i)
        pro_cds_dict_cowpox[name] = cowpox_cds_seq[3*i:3*(i+1)]
    pro_cds_dict_smallpox = {}
    smallpox_pro_seq = str(smallpox_pro.seq)
    smallpox_cds_seq = str(smallpox_cds.seq)
    for i in range(len(smallpox_pro_seq)):
        name = smallpox_pro_seq[i]+str(i)
        pro_cds_dict_smallpox[name] = smallpox_cds_seq[3*i:3*(i+1)]
    pro_cds_dict_monkeypox = {}
    monkeypox_pro_seq = str(monkeypox_pro.seq)
    monkeypox_cds_seq = str(monkeypox_cds.seq)
    for i in range(len(monkeypox_pro_seq)):
        name = monkeypox_pro_seq[i]+str(i)
        pro_cds_dict_monkeypox[name] = monkeypox_cds_seq[3*i:3*(i+1)]
    pro_cds_dict_skunkpox = {}
    skunkpox_pro_seq = str(skunkpox_pro.seq)
    skunkpox_cds_seq = str(skunkpox_cds.seq)
    for i in range(len(skunkpox_pro_seq)):
        name = skunkpox_pro_seq[i]+str(i)
        pro_cds_dict_skunkpox[name] = skunkpox_cds_seq[3*i:3*(i+1)]
    i = 0
    for seq_record in protein_seq:
        if seq_record.description.count("NC_006998")> 0:
            vaccinia_pro_align = seq_record
        elif seq_record.description.count("Cowpox")> 0:
            cowpox_pro_align = seq_record
        elif seq_record.description.count("Variola")> 0:
            smallpox_pro_align = seq_record
        elif seq_record.description.count("Monkey") > 0:
            monkeypox_pro_align = seq_record
        elif seq_record.description.count("Skunkpox") > 0:
            skunkpox_pro_align = seq_record
        i += 1
    i = 0
    vaccinia_cds_align = ''
    for j in range(len(vaccinia_pro_align.seq)):
        if str(vaccinia_pro_align.seq)[j] != '-':
            name = str(vaccinia_pro_align.seq)[j]+str(i)
            vaccinia_cds_align += pro_cds_dict_vaccinia[name]
            i += 1
        elif str(vaccinia_pro_align.seq)[j] == '-':
            vaccinia_cds_align += '---'
    i = 0
    cowpox_cds_align = ''
    for j in range(len(cowpox_pro_align.seq)):
        if str(cowpox_pro_align.seq)[j] != '-':
            name = str(cowpox_pro_align.seq)[j]+str(i)
            cowpox_cds_align += pro_cds_dict_cowpox[name]
            i += 1
        elif str(cowpox_pro_align.seq)[j] == '-':
            cowpox_cds_align += '---'
    i = 0
    smallpox_cds_align = ''
    for j in range(len(smallpox_pro_align.seq)):
        if str(smallpox_pro_align.seq)[j] != '-':
            name = str(smallpox_pro_align.seq)[j]+str(i)
            smallpox_cds_align += pro_cds_dict_smallpox[name]
            i += 1
        elif str(smallpox_pro_align.seq)[j] == '-':
            smallpox_cds_align += '---'
    i = 0
    monkeypox_cds_align = ''
    for j in range(len(monkeypox_pro_align.seq)):
        if str(monkeypox_pro_align.seq)[j] != '-':
            name = str(monkeypox_pro_align.seq)[j]+str(i)
            monkeypox_cds_align += pro_cds_dict_monkeypox[name]
            i += 1
        elif str(monkeypox_pro_align.seq)[j] == '-':
            monkeypox_cds_align += '---'
    i = 0
    skunkpox_cds_align = ''
    for j in range(len(skunkpox_pro_align.seq)):
        if str(skunkpox_pro_align.seq)[j] != '-':
            name = str(skunkpox_pro_align.seq)[j]+str(i)
            skunkpox_cds_align += pro_cds_dict_skunkpox[name]
            i += 1
        elif str(skunkpox_pro_align.seq)[j] == '-':
            skunkpox_cds_align += '---'
    fw_temp = open("temp.seq",'w')
    fw_temp.write("  5  "+str(len(vaccinia_cds_align)))
    fw_temp.write("\n")
    fw_temp.write("vaccinia")
    fw_temp.write("\n")
    fw_temp.write(vaccinia_cds_align)
    fw_temp.write("\n")
    fw_temp.write("cowpox")
    fw_temp.write("\n")
    fw_temp.write(cowpox_cds_align)
    fw_temp.write("\n")
    fw_temp.write("smallpox")
    fw_temp.write("\n")
    fw_temp.write(smallpox_cds_align)
    fw_temp.write("\n")
    fw_temp.write("monkeypox")
    fw_temp.write("\n")
    fw_temp.write(monkeypox_cds_align)
    fw_temp.write("\n")
    fw_temp.write("skunkpox")
    fw_temp.write("\n")
    fw_temp.write(skunkpox_cds_align)
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
fr = open("vaccinia_tree.txt")
fw = open("vaccinia_tree_process.txt",'w')
for lines in fr:
    fw.write(lines.replace('_',''))
fr.close()
fw.close()
