from sys import path, argv
import os
#path.append(f'{os.path.dirname(argv[0])}/../src/')
import SPlib as library
from sklearn.decomposition import PCA 
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB import NeighborSearch
#from Bio.Data.SCOPData import protein_letters_3to1 #old version
from Bio.Data.PDBData import protein_letters_3to1_extended as protein_letters_3to1
import gzip
import copy

#DSSP_DICT = {'H': 1, 'G': 1, 'I': 1, 'E': 2, 'B': 2}
DSSP_DICT = {'H': 1, 'E': 2}
       
def rdDSSP(fn, seq):
    data = []
    data0 = []
    mask = []
    active = False
    i0 = 0
    with open(fn) as f:
        for i, line in enumerate(f):
            ss = line.split()
            if ss[0] == '#' and ss[1] == 'RESIDUE': active=True; continue
            if not active: continue
            if line[13]=='!': continue
            while (line[13]!=seq[i0]) and line[13]!='X': #fix chain breaks here
                data.append(0)
                data0.append('C')
                mask.append(1)
                i0+=1
            data.append(DSSP_DICT.get(line[16], 0)) 
            data0.append(line[16]) 
            i0+=1
            if data[-1] == 2 and line[21]=='S': #split kinked strand
                mask.append(0)
            else:
                mask.append(1)
    while len(data)<len(seq): #zero resolution residues at terminus are not fixed by chain break above eg d1ml9a_
        data.append(0)
        data0.append('C')
        mask.append(1)
    return data, data0, mask


def rd(fn):
    with open(fn) as f:
        for line in f:
            ss = line.split()
            yield ss[0]

def get_axis(coords):
    center = coords.mean(axis=0)
    centered_coords = coords - center

    pca = PCA(n_components=3, svd_solver='full')
    pca.fit(centered_coords)

    principal_axis = pca.components_[0,:]
    return center, principal_axis

def get_proj(p, v, t): 
    p_translated = p - t 
    dot_product = np.dot(p_translated, v)
    magnitude = np.linalg.norm(v)
    proj = dot_product / magnitude * v / magnitude  #reuse magnitude
    return proj + t 

def write_file(fn, seq, coord, dssp, ideal_coord, starts, mids, ends):
    # This is inefficient storage but convenient for now (no need to duplicate coordinates)
    with open(fn, 'w') as f:
        f.write("%d %d\n"%(len(seq), len(starts))) #nres, nseg
        for aa, c, d, i in zip(seq, coord, dssp, ideal_coord):
            f.write("%c %.3f %.3f %.3f %d %.3f %.3f %.3f\n"%(aa, c[0], c[1], c[2], d, i[0], i[1], i[2]))
        for s, m, e in zip(starts, mids, ends):
            f.write("%d %d %d\n"%(s, m, e))

def write_empty(fn):
    with open(fn, 'w') as f:
        f.write("%d %d\n"%(0, 0)); return #nres, nseg

def check_nonlocal(ref, atom_list):
    count = 0 
    ref_id = ref.get_parent().id[1]
    for atom in atom_list:
        atom_id = atom.get_parent().id[1]
        if abs(ref_id-atom_id)>8: count+=1
    return count

def smooth_contacts(contact_number, cutoff=4.5, cutoff2=0.7):
    window_size = min(25, contact_number.shape[0])
    window = np.ones(window_size, dtype=np.float32)/window_size
    
    # Flag residues if they are part of a forward/backward span with low contact number 
    # Padding wont influence outcome
    smoothed_contact = np.convolve(contact_number, window, 'valid')
    forward = np.pad(smoothed_contact, pad_width=(0,window_size-1), mode='constant', constant_values=[10])
    backward = np.pad(smoothed_contact, pad_width=(window_size-1,0), mode='constant', constant_values=[10])

    # Ugly hack to handle short cases with overlapping padding
    mask = np.zeros_like(smoothed_contact, dtype=np.int32)
    mask1 = np.pad(mask, pad_width=(0,window_size-1), mode='constant', constant_values=[1])
    mask2 = np.pad(mask, pad_width=(window_size-1,0), mode='constant', constant_values=[1])
    mask = np.logical_not(mask1 & mask2).astype(np.float32)
    forward = forward * mask
    backward = backward * mask

    discard = ((forward<cutoff) | (backward<cutoff))

    # Trim residues if they are in a region enriched with flagged residues
    #window_size = min(51, contact_number.shape[0])
    window_size = 51
    discard = np.pad(discard, pad_width=[window_size//2]*2, mode='reflect')
    #Check behaviour of reflect if length is less than window_size//2
    discard = np.convolve(discard, np.ones(window_size, dtype=np.float32)/window_size, 'valid')
    discard = discard>cutoff2

    keep = np.logical_not(discard)
    return keep

def smooth_contacts_afdb_clust(contact_number, cutoff=4.5, cutoff2=0.7):
    window_size = min(25, contact_number.shape[0])
    window = np.ones(window_size, dtype=np.float32)/window_size
    
    smoothed_contact = np.convolve(contact_number, window, 'valid')
    #intentionally above threshold since it shouldn't impact discard decision
    forward = np.pad(smoothed_contact, pad_width=(0,window_size-1), mode='constant', constant_values=[10])
    backward = np.pad(smoothed_contact, pad_width=(window_size-1,0), mode='constant', constant_values=[10])

    discard = (forward<cutoff) | (backward<cutoff)

    #window_size = min(51, contact_number.shape[0])
    window_size = 51
    discard = np.pad(discard, pad_width=[window_size//2]*2, mode='reflect')
    discard = np.convolve(discard, np.ones(window_size, dtype=np.float32)/window_size, 'valid')
    discard = discard>cutoff2

    keep = np.logical_not(discard)
    return keep

def read_structure(fn, trim=False, parser=None):
    if fn.split('.')[-1] == 'gz': open_fn = gzip.open
    else: open_fn = open
    if not parser: p = PDBParser()
    else: p = parser
    with open_fn(fn, 'rt') as g:
        s = p.get_structure("", g)[0]

        residues = list(s.get_residues())
        
        seq = "".join([protein_letters_3to1.get(res.resname, 'X') for res in residues if 'CA' in res])
        nr = len(seq)
        coord = [list(res['CA'].coord) for res in residues if 'CA' in res]
        
        if trim:
            atoms = [a for a in s.get_atoms() if a.name == 'CA']
            n = NeighborSearch(atoms)

            contact_number = []
            for atom in atoms:
                contact_number.append(check_nonlocal(atom, n.search(atom.coord, 12)))

            trim_mask = smooth_contacts(np.array(contact_number))
            return nr, seq, np.array(coord), trim_mask
        else:
            return nr, seq, np.array(coord), np.ones(nr)

def process_structure(name, dssdir, sdir, odir, structure_suffix='pdb', af2model=False, trim=False, parser=None):
    fn = f'{sdir}/{name}.{structure_suffix}.gz'
    if not os.path.isfile(fn): fn = f'{sdir}/{name}.{structure_suffix}'
    nr, seq, coord, trim_mask = read_structure(fn, trim=trim, parser=parser)
    
    # Extract pre-computed DSSP labels
    try:
        dssp, dssp0, mask = rdDSSP(dssdir+'/'+name+'.dssp', seq)
    except FileNotFoundError:
        dssp = [0 for _ in range(nr)]
        mask = np.ones_like(dssp)
    dssp = np.array(dssp)
    mask = np.array(mask)
    
    if trim:
        nr = np.sum(trim_mask)
        #print(nr, trim_mask.shape[0])
        #print(fn)
        #print(seq)
        seq = "".join(np.array(list(seq))[trim_mask])
        coord = coord[trim_mask]
        dssp = dssp[trim_mask]
        mask = mask[trim_mask]
        #print(coord.shape)
        if len(seq)==0: 
            write_empty(f"{odir}/{name}.ideal")
            return

    if not af2model:
        # Compute rough SS labels from local CA distances
        tpl_dummy = library.Protein2(nr, seq, coord) #Dummy used since segs not reset
        tpl_dummy.calSS()
        ssec1 = tpl_dummy.getssec()

        # Merge DSSP and rough SS assignments
        dssp = np.max([dssp, ssec1], axis=0)
    
    #Keep only residues with CA
    dssp = dssp * np.array(mask).astype(np.int32)

    # Construct SPfast protein object    
    tpl = library.Protein2(nr, seq, coord)

    # Define SS segments
    tpl.calSS2(dssp) #This may merge segments not separated by a labeled loop

    starts = tpl.getsegstart()
    mids = tpl.getsegmid()
    ends = tpl.getsegend()
    qss = tpl.getssec()
    #print(fn)
    #print(ends)
    coord0 = copy.deepcopy(coord)

    #coord0 contains copy of original coordinates
    #coord used to update representative pseudoatoms

    # Iterate over valid segments
    for i in range(len(starts)):
        qfrag = []
        idxs1 = (starts[i], mids[i], ends[i]) # Segment representative indices
        coords1 = coord[idxs1[0]:(idxs1[2]+1)] # Atom coordinates in segment

        ## NOTE: This characterization may not be ideal for short segments where the PC1 can angle across the segment

        # Project coordinates to first principal component as idealized representative pseudoatoms
        center, principal_axis = get_axis(coords1) 
        qfrag.append(get_proj(coords1[0], principal_axis, center))
        qfrag.append(center)
        qfrag.append(get_proj(coords1[-1], principal_axis, center))
        coord[idxs1[0]] = qfrag[0]
        coord[idxs1[1]] = qfrag[1]
        coord[idxs1[2]] = qfrag[2]

    #write_file(f"ideal/{name1}.ideal", seq, coord0, dssp, coord) #used this - should be the same as below
    write_file(f"{odir}/{name}.ideal", seq, coord0, qss, coord, starts, mids, ends)

if __name__ == '__main__':
    import os
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('input_list', help='list of PDB ids')
    parser.add_argument('--dssdir', default='DSSP/', help='directory with DSSP files')
    parser.add_argument('--sdir', default='structures/', help='directory with gzipped pdbs')
    parser.add_argument('--odir', default='ideal/', help='output directory')
    parser.add_argument('--structure_suffix', default='pdb')
    parser.add_argument('--af2model', action='store_true', help='Flag for AF2 models to avoid loops being treated as sheets')
    parser.add_argument('--trim', action='store_true', help='Flag to trim hyperexposed residues likely to be disordered')
    args = parser.parse_args()

    if 'cif' in args.structure_suffix: p = MMCIFParser()
    else: p = PDBParser()
    for target in rd(args.input_list):
        if os.path.isfile(f"{args.odir}/{target}.ideal"): continue
        #if not os.path.isfile(f"{args.sdir}{target}.{args.structure_suffix}.gz"): continue
        process_structure(target, args.dssdir, args.sdir, args.odir, \
            structure_suffix=args.structure_suffix, trim=args.trim, af2model=args.af2model, parser=p)


