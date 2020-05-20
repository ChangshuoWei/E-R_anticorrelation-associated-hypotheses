gencode = {
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
      'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
      'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}


fr = open("orf_seq_7310.fasta")
fw = open("orf_protein.fasta",'w')

for lines in fr :
    if lines.startswith(">"):
        seq_name = lines
    else :
        if len(lines.strip()) != "orf_protein_length" :
            continue
        aa_seq = ''
        nucl_seq = lines.strip().upper()
        for i in range(0,"orf_protein_length"-3,3):
            codon = nucl_seq[i:(i+3)]
            if codon.count("A")+codon.count("C")+codon.count("G")+codon.count("T") ==3:
                aa = gencode[codon]
            else:
                aa = 'X'
            aa_seq += aa
        fw.write(seq_name)
        fw.write(aa_seq)
        fw.write("\n")
fr.close()
fw.close()

fr = open("orf_protein.fasta")


list_row = []

for lines in fr :
    if lines.startswith(">"):
        continue
    else:
        list_row.append(lines.strip())
fr.close()
aa_length = len(list_row[0])
list_col = []
for i in range(aa_length):
    list_col.append("")
    for seq in list_row:
        list_col[i] += seq[i]

var_count = 0
for i in range(aa_length):
    check_col = list_col[i].replace("X","")
    seq = {}
    for amino in check_col:
        seq[amino] = check_col.count(amino)
    amino_list = []
    for keys in seq.keys():
        amino_list.append(float(seq[keys]))
    if len(amino_list) == 1:
        continue
    else:
        start = 0
        for i in range(len(amino_list)):
            var_count += ((amino_list[i]*(len(check_col)-start-amino_list[i])))/(len(check_col)*(len(check_col)-1)/2)
            start += amino_list[i]
print(str(var_count))
print(str(var_count/"orf_protein_length"))
