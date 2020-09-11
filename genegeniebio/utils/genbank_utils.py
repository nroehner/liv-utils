'''
genegeniebio-utils (c) GeneGenie Bioinformatics Ltd. 2020

All rights reserved.

@author:  neilswainston
'''
# pylint: disable=too-many-arguments
from Bio import Seq, SeqIO, SeqFeature
from Bio.Alphabet import IUPAC


class GenbankBuilder():
    '''Class to build a Genbank record.'''

    def __init__(self, sequence, locus, accession, desc='', circular=True):
        seq = Seq.Seq(sequence, alphabet=IUPAC.ambiguous_dna)
        self.__record = SeqIO.SeqRecord(id=accession, seq=seq, name=locus,
                                        description=desc)

        if circular:
            self.__record.annotations['topology'] = 'circular'

    def add_feature(self, start, end, forward, feat_id, typ='misc_feature'):
        '''Add feature.'''
        location = SeqFeature.FeatureLocation(start, end, 1 if forward else -1)
        feature = SeqFeature.SeqFeature(location, type=typ)
        feature.qualifiers['label'] = feat_id
        self.__record.features.append(feature)

    def get_record(self):
        '''Get record.'''
        return self.__record

    def write_record(self, outfile):
        '''Write record.'''
        with open(outfile, 'w') as fle:
            SeqIO.write(self.__record, fle, 'genbank')
