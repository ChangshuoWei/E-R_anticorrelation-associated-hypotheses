from codon import caculate_codon

seq = {}
fr = open('wuhan1.fa')

for lines in fr :
    if lines.startswith(">"):
        name = lines.split()[0]
    else:
        seq[name] = lines.strip()
fr.close()

fr = open('check_var_once_test.txt')
fw = open('non_syn_once_test.txt','w')
fw2 = open('non_syn_statistc_once_test.txt','w')
orf1a_non = 0
orf1a_syn = 0
orf1b_non = 0
orf1b_syn = 0
S_non = 0
S_syn = 0
orf3a_non = 0
orf3a_syn = 0
E_non = 0
E_syn = 0
M_non = 0
M_syn = 0
orf6_non = 0
orf6_syn = 0
orf7a_non = 0
orf7a_syn = 0
orf8_non = 0
orf8_syn = 0
N_non = 0
N_syn = 0
orf1ab_non = 0
orf1ab_syn = 0 
for lines in fr:
    if lines.startswith("pos"):
        ncov_name = lines.strip()
        fw.write(ncov_name)
        fw.write('\n')
    else :
        pos = int(float(lines.strip().split('\t')[0]))
        if pos > 266 and pos < 13468:
            judge = (pos-266)%3
            if judge == 0:
                refer_codon = seq[name][(pos-1):(pos+2)]
                alter_codon = lines.strip().split('\t')[2]+seq[name][pos:(pos+2)]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 1:
                refer_codon = seq[name][(pos-2):(pos+1)]
                alter_codon = seq[name][(pos-2)]+lines.strip().split('\t')[2]+seq[name][pos]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 2:
                refer_codon = seq[name][(pos-3):pos]
                alter_codon = seq[name][(pos-3):(pos-1)] + lines.strip().split('\t')[2]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            fw.write(str(pos))
            fw.write('\t')
            fw.write(refer_amino)
            fw.write('\t')
            fw.write(alter_amino)
            fw.write('\n')
            if refer_amino == alter_amino:
                orf1a_syn += 1
            else :
                orf1a_non += 1
        elif pos > 13468 and pos < 21555:
            judge = (pos-13468)%3
            if judge == 0:
                refer_codon = seq[name][(pos-1):(pos+2)]
                alter_codon = lines.strip().split('\t')[2]+seq[name][pos:(pos+2)]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 1:
                refer_codon = seq[name][(pos-2):(pos+1)]
                alter_codon = seq[name][(pos-2)]+lines.strip().split('\t')[2]+seq[name][pos]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 2:
                refer_codon = seq[name][(pos-3):pos]
                alter_codon = seq[name][(pos-3):(pos-1)] + lines.strip().split('\t')[2]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            fw.write(str(pos))
            fw.write('\t')
            fw.write(refer_amino)
            fw.write('\t')
            fw.write(alter_amino)
            fw.write('\n')
            if refer_amino == alter_amino:
                orf1b_syn += 1
            else:
                orf1b_non += 1
        elif pos > 21563 and pos < 25384:
            judge = (pos-21563)%3
            if judge == 0:
                refer_codon = seq[name][(pos-1):(pos+2)]
                alter_codon = lines.strip().split('\t')[2]+seq[name][pos:(pos+2)]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 1:
                refer_codon = seq[name][(pos-2):(pos+1)]
                alter_codon = seq[name][(pos-2)]+lines.strip().split('\t')[2]+seq[name][pos]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 2:
                refer_codon = seq[name][(pos-3):pos]
                alter_codon = seq[name][(pos-3):(pos-1)] + lines.strip().split('\t')[2]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            fw.write(str(pos))
            fw.write('\t')
            fw.write(refer_amino)
            fw.write('\t')
            fw.write(alter_amino)
            fw.write('\n')
            if refer_amino == alter_amino:
                S_syn += 1
            else:
                S_non += 1
        elif pos > 25393 and pos < 26220:
            judge = (pos-25393)%3
            if judge == 0:
                refer_codon = seq[name][(pos-1):(pos+2)]
                alter_codon = lines.strip().split('\t')[2]+seq[name][pos:(pos+2)]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 1:
                refer_codon = seq[name][(pos-2):(pos+1)]
                alter_codon = seq[name][(pos-2)]+lines.strip().split('\t')[2]+seq[name][pos]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 2:
                refer_codon = seq[name][(pos-3):pos]
                alter_codon = seq[name][(pos-3):(pos-1)] + lines.strip().split('\t')[2]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            fw.write(str(pos))
            fw.write('\t')
            fw.write(refer_amino)
            fw.write('\t')
            fw.write(alter_amino)
            fw.write('\n')
            if refer_amino == alter_amino:
                orf3a_syn += 1
            else:
                orf3a_non += 1
        elif pos > 26245 and pos < 26472:
            judge = (pos-26245)%3
            if judge == 0:
                refer_codon = seq[name][(pos-1):(pos+2)]
                alter_codon = lines.strip().split('\t')[2]+seq[name][pos:(pos+2)]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 1:
                refer_codon = seq[name][(pos-2):(pos+1)]
                alter_codon = seq[name][(pos-2)]+lines.strip().split('\t')[2]+seq[name][pos]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 2:
                refer_codon = seq[name][(pos-3):pos]
                alter_codon = seq[name][(pos-3):(pos-1)] + lines.strip().split('\t')[2]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            fw.write(str(pos))
            fw.write('\t')
            fw.write(refer_amino)
            fw.write('\t')
            fw.write(alter_amino)
            fw.write('\n')
            if refer_amino == alter_amino:
                E_syn += 1
            else:
                E_non += 1
        elif pos > 26523 and pos < 27191:
            judge = (pos-26523)%3
            if judge == 0:
                refer_codon = seq[name][(pos-1):(pos+2)]
                alter_codon = lines.strip().split('\t')[2]+seq[name][pos:(pos+2)]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 1:
                refer_codon = seq[name][(pos-2):(pos+1)]
                alter_codon = seq[name][(pos-2)]+lines.strip().split('\t')[2]+seq[name][pos]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 2:
                refer_codon = seq[name][(pos-3):pos]
                alter_codon = seq[name][(pos-3):(pos-1)] + lines.strip().split('\t')[2]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            fw.write(str(pos))
            fw.write('\t')
            fw.write(refer_amino)
            fw.write('\t')
            fw.write(alter_amino)
            fw.write('\n')
            if refer_amino == alter_amino:
                M_syn += 1
            else:
                M_non += 1
        elif pos > 27202 and pos < 27387:
            judge = (pos-27202)%3
            if judge == 0:
                refer_codon = seq[name][(pos-1):(pos+2)]
                alter_codon = lines.strip().split('\t')[2]+seq[name][pos:(pos+2)]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 1:
                refer_codon = seq[name][(pos-2):(pos+1)]
                alter_codon = seq[name][(pos-2)]+lines.strip().split('\t')[2]+seq[name][pos]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 2:
                refer_codon = seq[name][(pos-3):pos]
                alter_codon = seq[name][(pos-3):(pos-1)] + lines.strip().split('\t')[2]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            fw.write(str(pos))
            fw.write('\t')
            fw.write(refer_amino)
            fw.write('\t')
            fw.write(alter_amino)
            fw.write('\n')
            if refer_amino == alter_amino:
                orf6_syn += 1
            else:
                orf6_non += 1
        elif pos > 27394 and pos < 27759:
            judge = (pos-27394)%3
            if judge == 0:
                refer_codon = seq[name][(pos-1):(pos+2)]
                alter_codon = lines.strip().split('\t')[2]+seq[name][pos:(pos+2)]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 1:
                refer_codon = seq[name][(pos-2):(pos+1)]
                alter_codon = seq[name][(pos-2)]+lines.strip().split('\t')[2]+seq[name][pos]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 2:
                refer_codon = seq[name][(pos-3):pos]
                alter_codon = seq[name][(pos-3):(pos-1)] + lines.strip().split('\t')[2]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            fw.write(str(pos))
            fw.write('\t')
            fw.write(refer_amino)
            fw.write('\t')
            fw.write(alter_amino)
            fw.write('\n')
            if refer_amino == alter_amino:
                orf7a_syn += 1
            else:
                orf7a_non += 1
        elif pos > 27894 and pos < 28259:
            judge = (pos-27894)%3
            if judge == 0:
                refer_codon = seq[name][(pos-1):(pos+2)]
                alter_codon = lines.strip().split('\t')[2]+seq[name][pos:(pos+2)]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 1:
                refer_codon = seq[name][(pos-2):(pos+1)]
                alter_codon = seq[name][(pos-2)]+lines.strip().split('\t')[2]+seq[name][pos]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 2:
                refer_codon = seq[name][(pos-3):pos]
                alter_codon = seq[name][(pos-3):(pos-1)] + lines.strip().split('\t')[2]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            fw.write(str(pos))
            fw.write('\t')
            fw.write(refer_amino)
            fw.write('\t')
            fw.write(alter_amino)
            fw.write('\n')
            if refer_amino == alter_amino:
                orf8_syn += 1
            else:
                orf8_non += 1
        elif pos > 28274 and pos < 29533:
            judge = (pos-28274)%3
            if judge == 0:
                refer_codon = seq[name][(pos-1):(pos+2)]
                alter_codon = lines.strip().split('\t')[2]+seq[name][pos:(pos+2)]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 1:
                refer_codon = seq[name][(pos-2):(pos+1)]
                alter_codon = seq[name][(pos-2)]+lines.strip().split('\t')[2]+seq[name][pos]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            elif judge == 2:
                refer_codon = seq[name][(pos-3):pos]
                alter_codon = seq[name][(pos-3):(pos-1)] + lines.strip().split('\t')[2]
                refer_amino = caculate_codon(refer_codon)
                alter_amino = caculate_codon(alter_codon)
            fw.write(str(pos))
            fw.write('\t')
            fw.write(refer_amino)
            fw.write('\t')
            fw.write(alter_amino)
            fw.write('\n')
            if refer_amino == alter_amino:
                N_syn += 1
            else:
                N_non += 1
orf1ab_non = orf1a_non + orf1b_non
orf1ab_syn = orf1a_syn + orf1b_syn
fw2.write('')
fw2.write('\t')
fw2.write('NS')
fw2.write('\t')
fw2.write('S')
fw2.write('\n')
fw2.write('orf1ab')
fw2.write('\t')
fw2.write(str(orf1ab_non))
fw2.write('\t')
fw2.write(str(orf1ab_syn))
fw2.write('\n')
fw2.write('S')
fw2.write('\t')
fw2.write(str(S_non))
fw2.write('\t')
fw2.write(str(S_syn))
fw2.write('\n')
fw2.write('orf3a')
fw2.write('\t')
fw2.write(str(orf3a_non))
fw2.write('\t')
fw2.write(str(orf3a_syn))
fw2.write('\n')
fw2.write('E')
fw2.write('\t')
fw2.write(str(E_non))
fw2.write('\t')
fw2.write(str(E_syn))
fw2.write('\n')
fw2.write('M')
fw2.write('\t')
fw2.write(str(M_non))
fw2.write('\t')
fw2.write(str(M_syn))
fw2.write('\n')
fw2.write('orf6')
fw2.write('\t')
fw2.write(str(orf6_non))
fw2.write('\t')
fw2.write(str(orf6_syn))
fw2.write('\n')
fw2.write('orf7a')
fw2.write('\t')
fw2.write(str(orf7a_non))
fw2.write('\t')
fw2.write(str(orf7a_syn))
fw2.write('\n')
fw2.write('orf8')
fw2.write('\t')
fw2.write(str(orf8_non))
fw2.write('\t')
fw2.write(str(orf8_syn))
fw2.write('\n')
fw2.write('N')
fw2.write('\t')
fw2.write(str(N_non))
fw2.write('\t')
fw2.write(str(N_syn))
fw2.write('\n')
fr.close()
fw.close()
fw2.close()
