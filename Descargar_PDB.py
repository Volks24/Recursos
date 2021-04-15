### Descargar PDB's
import Bio
PDB_All =['5T1A','6MEO','5WB2','5WB1','6LFL]

from Bio.PDB import PDBList
'''Selecting structures from PDB'''
pdbl = PDBList()

for file in PDB_All:
    pdbl.retrieve_pdb_file(file,pdir='PDB')
