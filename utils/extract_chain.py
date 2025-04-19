#!/bin/env python

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB import PDBIO
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.PDBIO import Select
from sys import argv

class ChainSelect(Select):
    def __init__(self, chain):
        self.chain = chain

    def accept_chain(self, chain):
        return chain.id == self.chain

fpath = argv[1].split('.')
if fpath[-1] == 'pdb':
    p = PDBParser()
    io = PDBIO()
    bn = argv[1].split('.pdb')[0]
    extension = 'pdb'
if fpath[-1] == 'cif':
    p = MMCIFParser()
    io = MMCIFIO()
    bn = argv[1].split('.cif')[0]
    extension = 'cif'
else:
    raise NotImplementedError

s = p.get_structure("", argv[1])

if len(argv)>2:
    chain_id = argv[2]
    chain_select = ChainSelect(chain_id)
    out_fn = f"{bn}_{chain_id}.{extension}"
else:
    chain_id = next(s.get_chains()).id
    chain_select = ChainSelect(chain_id)
    out_fn = argv[1]

io.set_structure(s)
io.save(out_fn, select=chain_select)

