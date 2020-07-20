'''
SBOL read example.
'''
from sbol import Document


def read(filename):
    '''Parses SBOL.'''
    doc = Document()
    doc.read(filename)

    comp_def = doc.componentDefinitions[0]

    # Get start / end:
    for annot in comp_def.sequenceAnnotations:
        location = annot.locations
        range = location.getRange()


if __name__ == '__main__':
    read('sbol.xml')
