from Bio import SeqIO
from Bio import SeqUtils
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
import pandas as pd

path = '/path/to/work/directory/'
infile = path + 'sequence.gb'
# record = SeqIO.parse(infile, "genbank")
record = SeqIO.read(open(infile, "r"), "genbank")
#rec = next(record)
recList = []
featureList = []
for feature in record.features:
    if feature.type == 'CDS':
        featureDic = feature.qualifiers
        if 'gene' in featureDic:
            gene = featureDic['gene']
        else:
            gene = featureDic['locus_tag']
        product = featureDic['product']
        start = feature.location.start.position
        end = feature.location.end.position
        seq = record.seq[start:end]
        if seq[:3] != 'ATG':
            seq = seq.reverse_complement()
        length = end - start
        translation = feature.qualifiers['translation'][0]
        mw = round(SeqUtils.molecular_weight(translation, seq_type='protein')/1000, 1)
        isoP = round(IP(translation).pi())
        featureList.append([gene, product, start, end, seq, length, translation, mw, isoP])
        df = pd.DataFrame(featureList)
        df.columns = ['Gene', 'Product', 'Start', 'End', 'DNA Sequence', 'DNA length', 'Protein Sequence',
                      'Molecular weight', 'Isoelectric point']
        df.to_excel(path + 'vacciniaGenome.xlsx')
