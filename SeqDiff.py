# This script picks different positions in 2 aligned sequences

from Bio import SeqIO
import pandas as pd

path = '/path/to/work/directory/'
infile = path + 'input.fasta'
seqList = []
records = SeqIO.parse(infile, 'fasta')
for record in records:
    seqName = record.description
    seq = list(record.seq._data.upper().decode())
    seq = [seqName] + seq
    seqList.append(seq)
df = pd.DataFrame(seqList)
for i in range(1, df.shape[1]):
    if df.iloc[1, i] != '-' and df.iloc[0, i] != df.iloc[1, i]:
        print('position, theoretical, obtained: ', i, df.iloc[1, i], df.iloc[0, i])
