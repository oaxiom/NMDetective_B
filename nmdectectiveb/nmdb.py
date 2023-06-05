
import sys, os, gzip
from shlex import split as shlexsplit
from collections import Counter
import matplotlib.pyplot as plot

'''

Reimplementation of the NMDetective_B algorithm in https://www.nature.com/articles/s41588-019-0517-5

Lindeboom et al., 2019 Nat Gen.

Adapted from Babarinde et al., 2021 NAR.

This version is generalised for GTFs, and standalone commandline

'''

class NMDB:
    def __init__(self, log=None):
        assert log, 'You must provide a log handle'
        self.log = log

    def gtf_parser(self, filename, gziped=False):
        '''
        Simple parser for GTF
        '''

        if filename.endswith('.gz'):
            gziped = True

        if gziped:
            gtfoh = gzip.open(filename, "rt")
        else:
            gtfoh = open(filename, "rt")

        for lin in gtfoh:
            if lin[0] == '#':
                continue

            l = lin.strip().split('\t')

            # Try to do some quick fails to speed up
            if l[2] == 'gene': # No useful data here
                continue

            tid = l[8].find('transcript_id')
            if 'exon_num' in l[8]:
                enu = l[8].find('exon_num')
                exon_num = l[8][enu:enu+17].split(' ')[1].strip(';').strip('"')
            else:
                exon_num = 0

            r = {'chrom': l[0],
                'left': int(l[3]),
                'right': int(l[4]),
                'feature': l[2],
                'strand': l[6],
                'gene_type': 'protein_coding',
                'transcript_id': l[8][tid:].split(' ')[1].strip(';').strip('"'),
                'exon_num': int(exon_num),
                }

            '''
            # This is very slow, and as I just need the transcript_id I hack it out above;

            # Split the attributes column
            attrs = {}
            for item in l[8].strip().split("; "):
                if item:
                    if '
                    item = item.strip()
                    ss = shlexsplit(item)
                    key = ss[0]
                    value = ss[1].strip('"')
                    attrs[key] = value
            '''

            yield r

        return


    def NMDetective_B_score(self, inlastexon, orflength, exonlength, within_50nt_of_lastEJ):
        '''

        Scoring decision tree

        '''
        if inlastexon:
            return 0.0, 'Last exon\n(Score=0.0)'

        if orflength < 150: # Distance from STOP to START
            return 0.12, 'Start-proximal\n(Score=0.12)'

        if exonlength >407: # exon length the STOP codon is in
            return 0.41, 'Long exon\n(Score=0.41)'

        if within_50nt_of_lastEJ: # EJC will block NMD if STOP within 50 nt.
            return 0.2, '50 nt rule\n(Score=0.20)'

        return 0.65, 'Trigger NMD\n(Score=0.65)'

    def score(self, gtf_filename):
        self.log.info('Starting score')

        scores = []
        cats = []
        tid_table = {}

        # STATS
        __skipped_no_start_stop = 0

        # First make bundles for each transcript
        transcript_bundles = {}
        exon_num = 0
        for idx, item in enumerate(self.gtf_parser(gtf_filename)):
            if idx > 0 and idx % 1e6 == 0:
                self.log.info(f'{idx:,} done')

            if item['gene_type'] != 'protein_coding':
                continue

            if item['feature'] not in set(['exon', 'start_codon', 'stop_codon']):
                continue

            if item['feature'] == 'exon':
                exon_num += 1

            if item['transcript_id'] not in transcript_bundles:
                transcript_bundles[item['transcript_id']] = []

            transcript_bundles[item['transcript_id']].append(item)

        self.log.info(f'Processed {idx:,} entries in the GTF')

        for transcript_id, transcript in transcript_bundles.items():
            exonStarts = []
            exonEnds   = []
            exonNums   = []
            # I need to convert the transcript bundle into:
            chrom = transcript[0]['chrom']
            strand = transcript[0]['strand']
            START = 0
            STOP = 0

            for entry in transcript:
                if entry['feature'] == 'start_codon': # This is not quite right as the codon is 3 bp. However, I think for these purposes it's close enough
                    START = entry['left']
                elif entry['feature'] == 'stop_codon':
                    STOP = entry['left']

                elif entry['feature'] == 'exon':
                    if strand == '+':
                        exonStarts.append(entry['left'])
                        exonEnds.append(entry['right'])
                    else:
                        exonStarts.append(entry['right'])
                        exonEnds.append(entry['left'])
                    exonNums.append(entry['exon_num'])

            exonStarts.sort() # seems they can be any order in the GFT. So put them in a specific order that I expect.
            exonEnds.sort()

            total_exons = len(exonStarts)

            # I need to work out the orflength
            orflength = 0 # currently wrong...
            for exons in zip(exonStarts, exonEnds):
                if strand == '+':
                    if START >= exons[0] and START <= exons[1] and STOP >= exons[0] and STOP <= exons[1]:
                        # START and STOP are in the same exon
                        orflength += STOP - START

                    elif START >= exons[0] and START <= exons[1]:
                        # START is inside this exon
                        orflength +=  exons[1] - START

                    elif STOP >= exons[0] and STOP <= exons[1]:
                        # STOP is inside this exon
                        orflength += STOP - exons[0]

                    else:
                        # check it's not a UTR
                        if exons[0] >= START and exons[1] <= STOP:
                            orflength += exons[1] - exons[0]

                else:
                    if START >= exons[1] and START <= exons[0] and STOP >= exons[1] and STOP <= exons[0]:
                        # START and STOP are in the same exon
                        orflength += START - STOP

                    elif START >= exons[1] and START <= exons[0]:
                        # START is inside this exon
                        orflength += START - exons[1]

                    elif STOP >= exons[1] and STOP <= exons[0]:
                        # STOP is inside this exon
                        orflength +=  exons[0] - STOP

                    else:
                        # check it's not a UTR
                        if exons[0] <= START and exons[1] >= STOP:
                            orflength += exons[0] - exons[1]

            if START == 0 or STOP == 0:
                #self.log.info(f'{transcript_id} has no start_codon or stop_codon key in the GTF, skipping')
                __skipped_no_start_stop += 1
                continue

            inlastexon = False
            within_50nt_of_lastEJ = False

            for e0, e1, exon_num in zip(exonStarts, exonEnds, exonNums):
                #print(transcript)
                #print(exonStarts, exonEnds)
                #print(STOP, exons, exon_num, total_exons)
                if STOP >= e0 and STOP <= e1:
                    if strand == '+':
                        d = STOP - e0
                        if exon_num == total_exons:
                            inlastexon = True
                        if exon_num == total_exons-1 and d < 50: # 50 bp of penultimate exon
                            within_50nt_of_lastEJ = True
                    else:
                        1/0 # Should be impossible to reach here;
                    exonlength = e1 - e0

                if STOP >= e1 and STOP <= e0:
                    if strand == '-':
                        d = e0 - STOP
                        if exon_num == 1:
                            inlastexon = True
                        if exon_num == 2 and d < 50: # 50 bp of penultimate exon
                            within_50nt_of_lastEJ = True
                    else:
                        1/0 # Should be impossible to reach here;
                    exonlength = e0 - e1

            #print(strand, exonlength)

            #print(inlastexon, orflength, exonlength, within_50nt_of_lastEJ)
            nmd_score, nmd_class = self.NMDetective_B_score(inlastexon, orflength, exonlength, within_50nt_of_lastEJ)
            #print(transcript_id, nmd_class, START, STOP)
            #print()

            cats.append(nmd_class)
            scores.append(nmd_score)
            tid_table[transcript_id] = (nmd_score, nmd_class.replace('\n', ' '))

        self.log.info(f'Found {len(transcript_bundles):,} transcripts')

        self.scores = scores
        self.cats = cats
        self.tid_table = tid_table

        self.log.info(f'Skipped as no identifiable START/STOP: {__skipped_no_start_stop:,} transcripts')

        return scores, cats

    def save_tid_table(self, filename):
        '''
        **Purpose**

            Save a table of transcript_id and NMD score

        '''

        oh = open(filename, 'wt')
        for tid in sorted(self.tid_table):
            oh.write(f'{tid}\t{self.tid_table[tid][0]}\t{self.tid_table[tid][1]}\n')

        oh.close()


    def plots(self, label):
        '''
        Output a few plots

        '''
        # Pie chart;

        # I want the categories to be in order:
        pie_summ = {
            'Last exon\n(Score=0.0)': 0.0,
            'Start-proximal\n(Score=0.12)': 0.0,
            '50 nt rule\n(Score=0.20)': 0.0,
            'Long exon\n(Score=0.41)': 0.0,
            'Trigger NMD\n(Score=0.65)': 0.0,
            }

        cols = [
            'lightgrey',
            'grey',
            'gold',
            'darkorange',
            'darkred',
            ]

        for c in self.cats:
            pie_summ[c] += 1

        fig = plot.figure(figsize=[3,3])
        ax = fig.add_subplot(111)
        ax.pie(pie_summ.values(), labels=pie_summ.keys(), autopct='%1.1f%%',
            textprops={'fontsize': 6},
            colors=cols,
            startangle=90,
            counterclock=False,)
        fig.savefig(f'{label}.pie.pdf')
        self.log.info('Drew Pie chart')


if __name__ == '__main__':
    import logging

    logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)-8s: %(message)s',
                    datefmt='%m-%d %H:%M')

    log = logging.getLogger('nmdetect')
    log.setLevel(logging.INFO)

    mpl_logger = logging.getLogger('matplotlib') # Bodge to silence the matplotlib logging
    mpl_logger.setLevel(logging.WARNING)

    n = NMDB(log=log)
    n.score('../test/hg38.gencode.v42.top100k.gtf.gz')
    n.save_tid_table('../test/hg38.gencode.v42.top100k.tid_table.tsv')
    n.plots('../test/hg38.gencode.v42.top100k')
