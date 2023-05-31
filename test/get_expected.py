'''

You need to download the file

hg38_NMDetectiveB_Lindeboom_et_al.v2.gtf.gz

from Figshare:

https://figshare.com/ndownloader/files/17357378

This script processes and extracts the NMMDB scores.

'''

import gzip
import matplotlib.pyplot as plot

oh = gzip.open('hg38_NMDetectiveB_Lindeboom_et_al.v2.gtf.gz', 'rt')

nmds = []
with_nmd_score = 0

for idx, line in enumerate(oh):
    if 'nmd_score' not in line:
        continue

    with_nmd_score += 1

    l = line.strip().split('\t')[8].split('; ')[-1].strip(';').lstrip('nmd_score ').split(',')
    #l = list(set(l))
    l = [float(i) for i in l]

    #print(l)

    nmds += l

oh.close()

print(f'Processed {idx:,} entries in the GTF')
print(f'Found {with_nmd_score:,} transcripts with a NMD score')
print(f'Found {len(nmds):,} scores')

# convert the scores to categories

pie_summ = {
    'Last exon\n(Score=0.0)': 0.0,
    'Start-proximal\n(Score=0.12)': 0.0,
    'Long exon\n(Score=0.41)': 0.0,
    '50 nt rule\n(Score=0.20)': 0.0,
    'Trigger NMD\n(Score=0.65)': 0.0,
    }

for s in nmds:
    if s == 0.0:
        pie_summ['Last exon\n(Score=0.0)'] += 1
    elif s == 0.12:
        pie_summ['Start-proximal\n(Score=0.12)'] += 1
    elif s == 0.41:
        pie_summ['Long exon\n(Score=0.41)'] += 1
    elif s == 0.20:
        pie_summ['50 nt rule\n(Score=0.20)'] += 1
    elif s == 0.65:
        pie_summ['Trigger NMD\n(Score=0.65)'] += 1
    else:
        print('No sum!', s)

fig = plot.figure(figsize=[2,2])
ax = fig.add_subplot(111)
ax.pie(pie_summ.values(), labels=pie_summ.keys(), autopct='%1.1f%%',
    textprops={'fontsize': 4})
fig.savefig(f'expected.hg38.knownGene.pdf')
print('Done')
