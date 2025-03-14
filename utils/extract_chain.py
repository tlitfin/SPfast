#!/bin/env python

from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB.PDBIO import Select
from sys import argv

class ChainSelect(Select):
    def __init__(self, chain):
        self.chain = chain

    def accept_chain(self, chain):
        return chain.id == self.chain

p = PDBParser()
io = PDBIO()

s = p.get_structure("", argv[1])

if len(argv)>2:
    chain_select = ChainSelect(argv[2])
    out_fn = f"{argv[1].split('.pdb')[0]}_{argv[2]}.pdb"
else:
    chain_select = ChainSelect(next(s.get_chains()).id)
    out_fn = argv[1]

io.set_structure(s)
io.save(out_fn, select=chain_select)

