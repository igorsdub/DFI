# IPython log file
"""
Fasta Convert
-------------

Suite of Tools for converting to and from the
fasta sequence.
"""
import pdbio
import glob


mapres = {'ALA': 'A',
          'CYS': 'C',
          'ASP': 'D',
          'GLU': 'E',
          'PHE': 'F',
          'GLY': 'G',
          'HIS': 'H',
          'ILE': 'I',
          'LYS': 'K',
          'LEU': 'L',
          'MET': 'M',
          'PRO': 'P',
          'ARG': 'R',
          'GLN': 'Q',
          'ASN': 'N',
          'SER': 'S',
          'THR': 'T',
          'TRP': 'W',
          'TYR': 'Y',
          'VAL': 'V'}


def get_seq(fname):
    """
    Get the structural sequence from a PDB file.

    Input
    -----
    fname: file
       name of pdb file
    """
    ATOMS = []
    pdbio.pdb_reader(fname, ATOMS, CAonly=True, noalc=True, Verbose=False)
    for a in ATOMS:
        yield mapres[a.res_name]


def format_seq(seq):
    """
    Format sequence so it is every 100 lines

    Input
    -----
    seq: str
       String of the protein's sequence.

    Output
    ------
    seq_str: str
       Formatted Sequence
    """
    i, seq_str = 0, ''
    for s in seq:
        seq_str = seq_str + s
        i += 1  # TODO use textwrap library instead
        if(i == 100):
            seq_str = seq_str + '\n'
            i = 0
    return seq_str


def fafsa_format(fname, outfileobj=None):
    """
    Converts a sequence in pdb format into
    fafsa format.

    Input
    -----
    fname: file
       Name of pdb file
    outfileobj: file object
       file object to write out to can be a filename
       or just stdout

    Output
    ------
    fafsafmt: str
       fafsa format of pdb file.
    """
    title = fname.split('_')[0]  # could lead to a bug
    seq_str = format_seq([x for x in get_seq(fname)])
    fafsafmt = """
>{}| | |len={}
{}""".format(title, len(seq_str), seq_str)
    if outfileobj:
        outfileobj.write(fafsafmt)
    else:
        print fafsafmt


def separate_fasta(fname):
    """
    Given a filename with multiple fasta sequences
    this function will separate them and write them
    to individual files labeled the title of the fasta
    sequence.fasta

    Input
    -----
    fname: file
       File containing multiple fasta sequences
    """
    outfile = None
    with open(fname, 'r') as infile:
        for line in infile:
            if line.startswith(">"):
                print "Begin"
                if outfile:
                    outfile.close()
                title = line[1:]
                title = str(title.strip(' '))
                print title
                outfile = open(title + '.fasta', 'w')
            outfile.write(line)
