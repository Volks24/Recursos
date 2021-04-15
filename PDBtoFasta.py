#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 16:54:59 2019
@author: sejmodha
"""
from Bio import Entrez,SeqIO

Entrez.email = "foo@bar.com"

query=r'1REV[All Fields] AND pdb[filter]'

handle=Entrez.esearch(db="protein", term=query)
records=Entrez.read(handle)
id_list=records['IdList']
#print(id_list)
handle.close()
for each_id in id_list:
    fasta=Entrez.efetch(db="protein", id=each_id, rettype="fasta")
    fasta_record=SeqIO.read(fasta, "fasta")
    print(f'>{fasta_record.id}|{fasta_record.description}\n{fasta_record.seq}')