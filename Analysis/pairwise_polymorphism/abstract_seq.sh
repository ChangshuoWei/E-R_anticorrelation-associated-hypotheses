makeblastdb -in orf_seq.fa -dbtype nucl -out orf_seq
blastn -query gisaid_0512_7310.fa  -db orf_seq -out N_seq_7310.fasta -outfmt "6 delim=@ qseqid qseq"
sed 's/hCoV/>hCoV/' orf_seq_7310.fasta  > orf_seq_7310_temp.fasta

sed 's/@/\n/' orf_seq_7310_temp.fasta > orf_seq_7310_temp2.fasta

rm -rf orf_seq_7310_temp.fasta orf_seq_7310.fasta
rename orf_seq_7310_temp2.fasta orf_seq_7310.fasta orf_seq_7310_temp2.fasta
