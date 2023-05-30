
import sys, os

'''

Reimplementation of the NMDetective_B algorithm in https://www.nature.com/articles/s41588-019-0517-5

Lindeboom et al., 2019 Nat Gen.

Adapted from Babarinde et al., 2021 NAR.

This version is generalised for GTFs, and standalone commandline

'''

class NMDB:
    def __init__(self):

    def gtf_parser(self, filename):
        '''
        Simple parser for GFT
        '''


    def NMDetective_B_score(inlastexon, orflength, exonlength, within_50nt_of_lastEJ):
        '''

        Scoring decision tree

        '''
        if inlastexon:
            return 0.0, 'Last exon'

        if orflength < 150: # Distance from STOP to START
            return 0.12, 'Start-proximal'

        if exonlength >407: # exon length the STOP codon is in
            return 0.41, 'Long exon'

        if within_50nt_of_lastEJ: # EJC will block NMD if STOP within 50 nt.
            return 0.2, '50 nt rule'

        return 0.65, 'Trigger NMD'

    def score(self, gtf_filename):

    res = {
        'MS Hit': [],
        'No MS hit': [],
        'Variant': [],
        'Matching': [],
        }

    data = {
        'Matching': matching_coding,
        'Variant': variant_coding,
        'MS Hit': has_peptide_hit,
        'No MS hit': no_peptide_hit,
        }

    __notfound = 0
    __nostop = 0

    for k in data:
        for transcript in data[k]:
            strand = transcript['strand']

            cds_key_to_use = None
            if 'cds_loc' in transcript and transcript['cds_loc']:
                cds_key_to_use = 'cds_loc'
            elif 'cds_genome_loc' in transcript and transcript['cds_genome_loc']:
                cds_key_to_use = 'cds_genome_loc'
            else:
                cds_key_to_use = 'cds_local_to_genome'

            orflength = abs(transcript[cds_key_to_use]['right'] - transcript[cds_key_to_use]['left'])

            if transcript[cds_key_to_use]['right'] == 0:
                # Ignore misannotated GENCODE CDSs
                continue

            if strand == '+':
                STOP = transcript[cds_key_to_use]['right']
            else:
                STOP = transcript[cds_key_to_use]['left']

            inlastexon = False
            within_50nt_of_lastEJ = False

            if STOP == 0:
                __nostop += 1
                continue

            #print(strand, list(zip(transcript['exonStarts'], transcript['exonEnds'])))
            for exon_num, exons in enumerate(zip(transcript['exonStarts'], transcript['exonEnds'])):
                if STOP >= exons[0] and STOP <= exons[1]:
                    if strand == '+':
                        d = exons[1] - STOP
                        if exon_num == len(transcript['exonEnds'])-1:
                            inlastexon = True
                            if d < 50:
                                within_50nt_of_lastEJ = True

                    else:
                        d = STOP - exons[0]
                        if exon_num == 0:
                            inlastexon = True
                            if d < 50:
                                within_50nt_of_lastEJ = True

                    exonlength = exons[1] - exons[0]

            nmd_score, nmd_class = NMDetective_B_score(inlastexon, orflength, exonlength, within_50nt_of_lastEJ)

            res[k].append(nmd_score)

        return res

    def plots(self, label):
        '''
        Output a few plots

        '''

        # Pie chart;
