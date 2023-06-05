'''

You need to download the file

hg38_NMDetectiveB_Lindeboom_et_al.v2.gtf.gz

from Figshare:

https://figshare.com/ndownloader/files/17357378

This script processes and extracts the NMMDB scores.

The gtf reported on the website has scores for every position, as if a STOP was introduced
at that point in the CDS. But to compare to our NMDetective we only need the score at
the actual reported STOP seen in the genuine transcript.

For + strand transcripts, that will be the last score.
But for - strand transcripts, that will be the first score

'''

import gzip
import matplotlib.pyplot as plot

oh = gzip.open('hg38_NMDetectiveB_Lindeboom_et_al.v2.gtf.gz', 'rt')

nmds = []
with_nmd_score = 0
tid_table = {}

transcript_bundle = {}

for idx, line in enumerate(oh):
    line = line.strip().split('\t')
    strand = line[6]
    gtf_decorators = line[8]

    tid = gtf_decorators.find('transcript_id')
    tid = gtf_decorators[tid:].split(' ')[1].strip(';').strip('"')

    if tid not in transcript_bundle:
        transcript_bundle[tid] = {
            'STOP': None,
            'NMDB_scores': [],
            'strand': strand,
            }

    if 'stop_codon' in line:
        if strand == '+':
            transcript_bundle[tid]['STOP'] = int(line[3])-2 # It doen't give a score for the actual STOP?
        else:
            transcript_bundle[tid]['STOP'] = int(line[4])+2

    if 'nmd_score' in gtf_decorators:
        l = gtf_decorators.split('; ')[-1].strip(';').lstrip('nmd_score ').split(',')
        l = [float(i) for i in l]
        transcript_bundle[tid]['NMDB_scores'].append({'l': int(line[3]), 'r': int(line[4]), 's': l})

#print(transcript_bundle)

# Making a PIE chart of all the nmd_scores in the file perfectly replicates the figure in fig 1d of the paper.
no_stop = 0
no_nmd_data = 0
valid_transcripts = 0
for tid in transcript_bundle:
    STOP = transcript_bundle[tid]['STOP']
    if not STOP: # Didn't fid a STOP codon
        no_stop += 1
        continue

    if not transcript_bundle[tid]['NMDB_scores']:
        no_nmd_data += 1
        continue

    score = None
    valid_transcripts += 1

    # find the exon the STOP is in;
    # I don't need to care about strand this way, just grap the STOP codon;
    for exon in transcript_bundle[tid]['NMDB_scores']:
        if STOP >= exon['l'] and STOP <= exon['r']:
            localpos = exon['r'] -STOP
            score = exon['s'][localpos]

    if score is not None:
        with_nmd_score += 1
        tid_table[tid] = score
        nmds.append(score)
    else:
        #print(score)
        #print(tid, transcript_bundle[tid])
        pass

oh.close()

print(f'Processed {idx:,} entries in the GTF')
print(f'Found {len(transcript_bundle):,} transcripts')
print(f'  Of which {no_stop:,} had no STOP')
print(f'  Of which {no_nmd_data:,} had no NMDB score')
print(f'Leaving {valid_transcripts:,} transcripts')
print(f'Assigned {with_nmd_score:,} transcripts with a NMD score')

print('Save TID table')
oh = open('expected.Lindeboom.tidtable.tsv', 'wt')
for tid in sorted(tid_table):
    oh.write(f'{tid}\t{tid_table[tid]}\n')
oh.close()

# convert the scores to categories

pie_summ = {
    'Last exon\n(Score=0.0)': 0.0,
    'Start-proximal\n(Score=0.12)': 0.0,
    '50 nt rule\n(Score=0.20)': 0.0,
    'Long exon\n(Score=0.41)': 0.0,
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

cols = [
    'lightgrey',
    'grey',
    'gold',
    'darkorange',
    'darkred',
    ]

fig = plot.figure(figsize=[3,3])
ax = fig.add_subplot(111)
ax.pie(pie_summ.values(), labels=pie_summ.keys(), autopct='%1.1f%%',
            textprops={'fontsize': 6},
            colors=cols,
            startangle=90,
            counterclock=False,)
fig.savefig(f'expected.Lindeboom.pdf')
fig.savefig(f'expected.Lindeboom.png')
print('Done')
