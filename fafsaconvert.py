# IPython log file

import pdbio
import glob 
fname = glob.glob('*.pdb')

        
mapres={'ALA':'A',
'CYS':'C',
'ASP':'D',
'GLU':'E',
'PHE':'F',
'GLY':'G',
'HIS':'H',
'ILE':'I',
'LYS':'K',
'LEU':'L',
'MET':'M',
'PRO':'P',
'ARG':'R',
'GLN':'Q',
'ASN':'N',
'SER':'S',
'THR':'T',
'TRP':'W',
'TYR':'Y',
'VAL':'V'}
    
def get_seq(fname):
    ATOMS=[]
    pdbio.pdb_reader(fname,ATOMS,CAonly=True,noalc=True,Verbose=False)
    for a in ATOMS:
        yield mapres[a.res_name]
        
def format_seq(seq):
    i,seq_str = 0,''
    for s in seq:
        seq_str = seq_str + s
        i+=1
        if(i == 100):
            seq_str = seq_str+'\n'
            i = 0 
    return seq_str

def fafsa_format(fname,outfileobj=None):
    title = fname.split('_')[0]
    seq_str = format_seq( [x for x in get_seq(fname) ] )
    fafsafmt="""
>{}| | |len={}
{}""".format(title,len(seq_str),seq_str)
    if outfileobj:
        outfileobj.write(fafsafmt)
    else:
        print fafsafmt 

