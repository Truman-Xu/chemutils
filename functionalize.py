## Frank Lab 2021
## University of Michigan

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.rdMolAlign import AlignMol
import copy
from itertools import product, repeat
from tqdm.contrib.concurrent import process_map
import os

def findHIds(mol, unique: bool = True):
    '''
    Find the hydrogens from the molecule
    If unique flag is set only symmetrically unique hydrogens will be returned
    Return a dict with the format {atom order : (H index, Neighbor Atom index)}
    '''
    ## find symmetrically equivalent hydrogens by canonical order
    ## if unique, no equivalent order will be given (unique rank for equivalent Hs)
    orders = list(Chem.CanonicalRankAtoms(mol, breakTies = not unique))
    ## find atom indices of hydrogens and the indices of their neighbor atoms
    H_neigh_ids = [(a.GetIdx(), a.GetNeighbors()[0].GetIdx()) 
                   for a in mol.GetAtoms() if a.GetSymbol() == 'H']
    ## return only one hydrogen per rank {order : (H index, Neighbor index)}
    return {orders[ids[0]] : ids for ids in H_neigh_ids}

def connectMols(lig, frag, id_ligH, id_ligNbr, id_fragH, id_fragNbr):
    """
    Connect two mols whose hydrogens are aligned into position
    Return a combined mol
    Function copied from here with modification
    https://github.com/molecularsets/moses/blob/master/moses/baselines/combinatorial.py
    """
    # align the fragament to the ligand 
    # fragment will move to the ligand to preserve docked poses
    _rmsd = AlignMol(frag, lig, 
                     atomMap=((id_fragH, id_ligNbr),(id_fragNbr, id_ligH)))
    # combine and make an editable mol object
    # ligand first, fragment second to preserve the atom idx of ligand
    combined = Chem.CombineMols(lig, frag)
    emol = Chem.EditableMol(combined)
    # form the bond
    emol.AddBond(id_ligNbr,
                 id_fragNbr + lig.GetNumAtoms(), 
                 # fragament atom idx increase by the num of atom in ligand
                 order=Chem.BondType.SINGLE)
    # remove the hydrogens from the ligand and the fragment
    emol.RemoveAtom(id_fragH + lig.GetNumAtoms())
    emol.RemoveAtom(id_ligH)
    # return the immutable mol object
    return emol.GetMol()

def functionalizeH(lig, frag, unique_H: bool = True):
    '''
    Functionalize the ligand with a fragment on all unique positions of hydrogen atoms
    If unqiue flag is set to False, all hydrogen positions on the ligand will be used
    for forming the bonds. 
    Return a list of functionalized ligand
    '''
    fused = []
    # get combination of positions of hydrogens from both lig and frag
    # only unique hydrogens position will be used on the fragament
    combos = product(findHIds(lig, unique_H).values(), findHIds(frag, True).values())
    for combo in combos:
        id_ligH, id_ligNbr = combo[0]
        id_fragH, id_fragNbr = combo[1]
        # create a new copy (deep copy) of the mols
        cur_lig = copy.deepcopy(lig)
        cur_frag = copy.deepcopy(frag)
        fused.append(connectMols(cur_lig, cur_frag, id_ligH, id_ligNbr, id_fragH, id_fragNbr)) 
    return fused
    
def _FuncHWrapper(args):
    return functionalizeH(*args[0], unique_H = args[1])

def enumFuncH(lig_list, frag_list, unique_H: bool = True):
    '''
    Enumerated all the functionalized molecules 
    from a list of ligand and a list of fragment
    Return a list of lists with format 
    [[lig1-frag1, lig1-frag2, ..., lig1-fragN],
     [lig2-frag1, lig2-frag2, ..., lig2-fragN],
     ...
     [ligN-frag1, ligN-frag2, ..., ligN-fragN]]
    '''
    args = list(zip(product(lig_list, frag_list), repeat(unique_H)))
    fused = process_map(
            _FuncHWrapper,
            args,
            chunksize = 5
        )
    return fused

def sdfEnumFuncH(lig_path, frag_path, out_path, unique_H: bool = True):
    '''
    Enumerated all the functionalized molecules 
    from .sdf ligand file and .sdf fragment file 
    Files can contain multiple ligands or fragments to be connected
    Save the fused molecules in sdf files with presevation of ligands 3D coordinates
    '''
    if os.path.exists(out_path): 
        ligs = [m for m in Chem.SDMolSupplier(lig_path, removeHs = False)]
        frags = [m for m in Chem.SDMolSupplier(frag_path, removeHs = False)]
        # generate a list of the functioanlized molecules
        results = enumFuncH(ligs, frags, unique_H)
        # create a list of tuples with ids of ligand and fragment pair 
        # for list of the functioanlized molecules
        ids = list(product(range(len(ligs)), range(len(frags))))
        # save the molecules by ligand_fragment indices
        for i, fused in zip(ids, results):
            file_path = os.path.join(out_path, 'combined_{}_{}.sdf'.format(*i))
            with open(file_path, 'w') as f:
                w = Chem.SDWriter(f)
                for mol in fused:
                    w.write(mol)
                w.close()
    else: 
        raise FileNotFoundError("Directory doesn't exist:", out_path)
        
if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser(description='Ligand-Fragment Functionalization')
    parser.add_argument('-l','--lig', type=str,
                        help='ligand .sdf file path')
    parser.add_argument('-f','--frag', type=str,
                        help='fragment .sdf file path')
    parser.add_argument('-o','--out', type=str,
                        help='directory path for file output')
    parser.add_argument('--unique', type=bool, default=True,
                        help='only form bond to symmetrically unique hydrogens. Default to True')
    args = parser.parse_args()
    
    if not os.path.exists(args.out):
        os.makedirs(args.out)
        print('Directory Made:',args.out)
    
    sdfEnumFuncH(args.lig, args.frag, args.out, args.unique)