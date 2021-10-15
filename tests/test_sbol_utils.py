'''
liv-utils (c) University of Liverpool. 2020

All rights reserved.

@author:  neilswainston
'''
import os
import tempfile
import unittest

from liv_utils import sbol_utils


class Test(unittest.TestCase):
    '''Test class for dna_utils.'''

    directory = os.path.join(
        os.path.dirname(
            os.path.realpath(__file__)
        ),
        'data'
    )

    def test(self):
        '''Tests round trip equality.'''
        dna1 = sbol_utils.read(os.path.join(self.directory, 'sbol.xml'))
        dna2 = round_trip(dna1)
        self.assertEqual(dna1, dna2)


def round_trip(dna):
    '''Writes / reads DNA object, via SBOL export / import.'''
    tmp = tempfile.NamedTemporaryFile()
    sbol_utils.write(dna, tmp.name)
    return sbol_utils.read(tmp.name)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
