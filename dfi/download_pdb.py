"""
Download PDB
============

Description
-----------
Download pdb is suite of tools for interfacing
with the protein data bank (PDB)
"""
def fetch_pdb(id,Verbose=False):
    """
    Download pdb file and write out to file id.pdb

    Input
    ------
    id: str
       4 letter pdb code

    Output
    ------
    id.pdb: file
       filename id.pdb 
    """
    import urllib
    url = 'http://www.rcsb.org/pdb/files/%s.pdb' % id
    with open(id+'.pdb','w') as outfile:
        outfile.write(urllib.urlopen(url).read())
    if(Verbose):
        print("Wrote out %s.pdb"%(id))
