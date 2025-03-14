from pymol import cmd, selector, stored
import pprint
import math
import copy
from sklearn.decomposition import PCA 
import numpy as np

protein_letters_1to3 = {
    "A": "Ala",
    "C": "Cys",
    "D": "Asp",
    "E": "Glu",
    "F": "Phe",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "K": "Lys",
    "L": "Leu",
    "M": "Met",
    "N": "Asn",
    "P": "Pro",
    "Q": "Gln",
    "R": "Arg",
    "S": "Ser",
    "T": "Thr",
    "V": "Val",
    "W": "Trp",
    "Y": "Tyr",
}

stored.protein_letters_3to1 = {value.upper(): key for key, value in protein_letters_1to3.items()}

stored.ss_dict = {"H": 1, "S": 2}

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

def idealize(tpl, coord):
    starts = tpl.getsegstart()
    mids = tpl.getsegmid()
    ends = tpl.getsegend()
    qss = tpl.getssec()
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
    return coord

def SPfast( sel1, sel2 ):
    """
    SPfast-based Alignment of two protein structures
    Adapted from cealign plugin

    Params:
    \@param sel1: (string) valid PyMol selection string of protein 1 to align	
    \@param sel2: (string) valid PyMol selection string of protein 2 to align

    Side-Effects:
    \@note: rotates and translates the objects (proteins) provided in the
    selections, sel1 and sel2, to represent the alignment.  Probably will
    also change their representation to more clearly show the aligned
    segments.
    """   
    
    # make the lists for holding coordinates
    # partial lists
    stored.sel1 = []
    stored.ss1 = []
    stored.seq1 = []

    stored.sel2 = []
    stored.ss2 = []
    stored.seq2 = []

    # full lists
    stored.mol1 = []
    stored.mol2 = []
 
    # now put the coordinates into a list
    # partials
 
    # -- REMOVE ALPHA CARBONS
    sel1 = sel1 + " and n. CA"
    sel2 = sel2 + " and n. CA"
    # -- REMOVE ALPHA CARBONS
 
    cmd.iterate_state(1, selector.process(sel1), "stored.sel1.append([x,y,z])")
    cmd.iterate_state(1, selector.process(sel2), "stored.sel2.append([x,y,z])")
    
    # full molecule
    mol1 = cmd.identify(sel1,1)[0][0]
    mol2 = cmd.identify(sel2,1)[0][0]
    
    # put all atoms from MOL1 & MOL2 into stored.mol1
    cmd.iterate_state(1, mol1, "stored.mol1.append([x,y,z])")
    cmd.iterate_state(1, mol2, "stored.mol2.append([x,y,z])")

    # Store SS
    cmd.iterate_state(1, sel1, "stored.ss1.append(stored.ss_dict.get(ss, 0))")
    cmd.iterate_state(1, sel2, "stored.ss2.append(stored.ss_dict.get(ss, 0))")

    # Store seq
    cmd.iterate_state(1, sel1, "stored.seq1.append(stored.protein_letters_3to1.get(resn,'X'))")
    cmd.iterate_state(1, sel2, "stored.seq2.append(stored.protein_letters_3to1.get(resn,'X'))")
    
    if ( len(stored.mol1) == 0 ):
        print("ERROR: Your first selection was empty.")
        return
    if ( len(stored.mol2) == 0 ):
        print("ERROR: Your second selection was empty.")
        return
    
    import SPlib
    import numpy as np

    query = SPlib.Protein2(len(stored.seq1), "".join(stored.seq1), stored.sel1)
    query.calSS()
    ssec1 = query.getssec()
    dssp1 = np.max([ssec1, stored.ss1], axis=0)
    query.calSS2(dssp1)
    query.add_ideal(len(stored.seq1),idealize(query, np.array(stored.sel1)))

    tpl = SPlib.Protein2(len(stored.seq2), "".join(stored.seq2), stored.sel2)
    tpl.calSS()
    ssec2 = tpl.getssec()
    dssp2 = np.max([ssec2, stored.ss2], axis=0)
    tpl.calSS2(dssp2)
    tpl.add_ideal(len(stored.seq2), idealize(tpl, np.array(stored.sel2)))

    salign = SPlib.Salign()
    print("SPscore:", salign.run_pairwise(query, tpl))

    a = np.zeros((4,3))
    SPlib.get_u(salign, a)
    b = np.zeros((4,4))
    b[:3,:3] = a[:3,:3]
    b[:3,3] = a[3,:3]

    print("Transformation Matrix:", b)

    cmd.transform_selection(mol1, b.flatten(), homogenous=0)
    
## Let PyMol know about the command
cmd.extend("SPfast", SPfast)
