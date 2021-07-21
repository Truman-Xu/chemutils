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

def setTetherIdx(mol):
    mol2 = Chem.rdmolops.RemoveHs(mol)
    idx_list = [str(i) for i in range(1, mol2.GetNumAtoms()+1)]
    # set tethering atom ids for both versions of rdock
    mol.SetProp('rxdock.tethered_atoms',', '.join(idx_list))
    mol.SetProp('rdock.tethered_atoms',', '.join(idx_list))

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
    mol_out = emol.GetMol()
    Chem.SanitizeMol(mol_out)
    return mol_out

def functionalizeH(lig, frag, unique_H: bool = True, tether: bool = False ):
    '''
    Functionalize the ligand with a fragment on all unique positions of hydrogen atoms
    If unqiue flag is set to False, all hydrogen positions on the ligand will be used
    for forming the bonds. 
    Return a list of functionalized ligand
    '''
    # tether atom index set on the original ligand mol will be preserved when copied
    if tether: setTetherIdx(lig)
    fused_list = []
    # get combination of positions of hydrogens from both lig and frag
    # only unique hydrogens position will be used on the fragament
    combos = product(findHIds(lig, unique_H).values(), findHIds(frag, True).values())
    for combo in combos:
        id_ligH, id_ligNbr = combo[0]
        id_fragH, id_fragNbr = combo[1]
        # create a new copy (deep copy) of the mols
        cur_lig = copy.deepcopy(lig)
        cur_frag = copy.deepcopy(frag)
        fused = connectMols(cur_lig, cur_frag, id_ligH, id_ligNbr, id_fragH, id_fragNbr)
        fused_list.append(fused) 
    return fused_list
    
def _FuncHWrapper(args):
    return functionalizeH(*args[0], unique_H = args[1], tether=args[2])

def combinatorialFunc(lig_list, frag_list, unique_H: bool = True, tether: bool = False):
    '''
    Enumerated all the functionalized molecules 
    from a list of ligand and a list of fragment
    Return a 3D list of lists with format 
    [[lig1-frag1, lig1-frag2, ..., lig1-fragN],
     [lig2-frag1, lig2-frag2, ..., lig2-fragN],
     ...
     [ligN-frag1, ligN-frag2, ..., ligN-fragN]]
     where lig_i-frag_j is a list of mols fused on various hydrogen positions
    '''
    args = list(zip(product(lig_list, frag_list), repeat(unique_H), repeat(tether)))
    fused = process_map(
            _FuncHWrapper,
            args,
            chunksize = 100
        )
    return fused

def _sdfFuncWrapper(args):
    '''
    Function wrapper for multiprocessing on functionalization and writing to file
    each thread processes one ligand and write to file
    '''
    lig_idx, lig, frags, path_template, unique_H, tether = args
    # tether atom index set on the original ligand mol will be preserved when copied 
    if tether: setTetherIdx(lig)
    # get combination of positions of hydrogens from both lig and frag
    # only unique hydrogens position will be used on the fragament
    for frag_idx, frag in enumerate(frags):
        combos = product(findHIds(lig, unique_H).values(), findHIds(frag, True).values())
        file_path = path_template.format(lig_idx, frag_idx)
        with open(file_path, 'w') as f:
            w = Chem.SDWriter(f)
            for combo in combos:
                id_ligH, id_ligNbr = combo[0]
                id_fragH, id_fragNbr = combo[1]
                # create a new copy (deep copy) of the mols
                cur_lig = copy.deepcopy(lig)
                cur_frag = copy.deepcopy(frag)
                w.write(connectMols(cur_lig, cur_frag, id_ligH, id_ligNbr, id_fragH, id_fragNbr)) 
            w.close()

def sdfCombinatorialFunc(lig_path, frag_path, out_path, unique_H: bool = True, tether: bool = False):
    '''
    Enumerated all the combinatorial of ligands and fragments on all H positions and functionalize the pairs
    Takes input from .sdf ligand file and .sdf fragment file 
    Files can contain multiple ligands or fragments to be connected
    Save the fused molecules in sdf files with presevation of ligands 3D coordinates
    Multithreaded: each thread processes one ligand and write to files
    '''
    if os.path.exists(out_path): 
        # Determined the file path and set file name template
        lig_name = os.path.split(lig_path)[-1].split('.')[0]
        out_file_template = lig_name+'_combined_{}_{}.sdf'
        out_path_template = os.path.join(out_path, out_file_template)
        
        ligs = [m for m in Chem.SDMolSupplier(lig_path, removeHs = False)]
        frags = [m for m in Chem.SDMolSupplier(frag_path, removeHs = False)]

        args = list((i, lig, frags, out_path_template, unique_H, tether) for i, lig in enumerate(ligs))
        fused = process_map(
            _sdfFuncWrapper,
            args,
            chunksize = 1
        )
        print(len(ligs)*len(frags),"ligand-fragment pair processed.")
    else: 
        raise FileNotFoundError("Directory doesn't exist:", out_path)
        
if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser(description=
                                     '''
                                     Ligand-Fragment Functionalization
                                     DETAILS TO BE FILLED HERE
                                     '''
                                    )
    parser.add_argument('-l','--lig', type=str,
                        help='ligand .sdf file path')
    parser.add_argument('-f','--frag', type=str,
                        help='fragment .sdf file path')
    parser.add_argument('-o','--outpath', type=str, default="./out",
                        help='Directory path for file output. Default to the ./out')
    parser.add_argument('--allH', action='store_true', default=False,
                        help='''
                        Attach to all Hydrogen position instead of symmetrically unique ones. \n
                        Without this flag, it only replace symmetrically unique hydrogens
                        ''')
    parser.add_argument('--tether', action='store_true', default=False,
                        help='''
                        Set the tether docking atom ids for the functionalized ligand.\n
                        If the flag is present, the atoms from the ligand will be fixed during docking, 
                        and only the fragment parts are allowed to move freely.\n
                        Default to no tethering (without this flag)
                        ''')
    args = parser.parse_args()
    
    if not os.path.exists(args.outpath):
        os.makedirs(args.outpath)
        print('Directory Made:',args.outpath)
        
    uniqueH = not args.allH
    sdfCombinatorialFunc(args.lig, args.frag, args.outpath, uniqueH, args.tether)