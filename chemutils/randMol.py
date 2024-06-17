import random
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondType
from rdkit.Chem.rdMolAlign import AlignMol
from copy import deepcopy
from .fragmentation import PrepMolFrags, MakeMolFragDict, DrawMolFragDict

def GetAtomDict(mol: Chem.Mol):
    atom_idx_dict = {}
    for a in mol.GetAtoms():
        if a.GetSymbol() == '*' or a.GetSymbol() == 'H':
            assert len(a.GetNeighbors()) == 1
            Nei_idx = a.GetNeighbors()[0].GetIdx()
            if Nei_idx in atom_idx_dict:
                atom_idx_dict[Nei_idx].append(a.GetIdx())
            else: 
                atom_idx_dict[Nei_idx] = [a.GetIdx()]
    return atom_idx_dict

def PrepAttachPts(mol: Chem.RWMol):
    atom_idx_dict = GetAtomDict(mol)
    heavy_idx = random.choice(list(atom_idx_dict.keys()))
    wild_idx = random.choice(atom_idx_dict[heavy_idx])
    return wild_idx, heavy_idx

def PrepRandFrag(mol_fragment_list):
    choice_fragment = random.choice(mol_fragment_list)
    atom_idx_dict = GetAtomDict(choice_fragment)
    heavy_idx = random.choice(list(atom_idx_dict.keys()))
    wild_idx = random.choice(atom_idx_dict[heavy_idx])
    return choice_fragment, wild_idx, heavy_idx

def _RandAssemFrags(temp_frag_list, is3D):
    
    init_mol = random.choice(temp_frag_list)
    temp_frag_list.remove(init_mol)
    ed_mol = Chem.RWMol(init_mol)
    
    while len(temp_frag_list)>0:
        subs_idx1, heavy_idx1 = PrepAttachPts(ed_mol)
        
        choice_fragment, subs_idx2, heavy_idx2 = \
        PrepRandFrag(temp_frag_list)
        
        if is3D:
        # align the fragament to the ed_mol 
        # fragment will move to the ed_mol
            _rmsd = AlignMol(choice_fragment, ed_mol, 
                atomMap=(
                    (subs_idx2, heavy_idx1),
                    (heavy_idx2, subs_idx1)
                    )
                )
        # Num of atoms before insert the choice fragment
        cur_num_atoms = ed_mol.GetNumAtoms()

        ed_mol.InsertMol(choice_fragment)
        ed_mol.AddBond(
            heavy_idx1, 
            # fragament atom idx increase by the num of atom in ed_mol
            # in the previous state
            heavy_idx2+cur_num_atoms, 
            BondType.SINGLE
            )
        # remove the hydrogens or wild card from the ed_mol and the fragment
        # starting from the higher idx atom
        ed_mol.RemoveAtom(subs_idx2+cur_num_atoms)
        ed_mol.RemoveAtom(subs_idx1)
        
        temp_frag_list.remove(choice_fragment)

    for a in ed_mol.GetAtoms():
        if a.GetSymbol() == '*':
            ed_mol.ReplaceAtom(a.GetIdx(), Chem.Atom(1))

    final_mol = ed_mol.GetMol()
    final_mol.UpdatePropertyCache()
    Chem.SanitizeMol(final_mol)
    if is3D:
        AllChem.UFFOptimizeMolecule(final_mol)
    final_mol = Chem.RemoveHs(final_mol)

    return final_mol

def RandAssemFrags(fragment_dict, is3D):
    temp_frag_list = []
    for smiles, sub_dict in fragment_dict.items():
        for i in range(sub_dict['n']):
            temp_frag_list.append(sub_dict['mol'])
    try:
        return _RandAssemFrags(temp_frag_list, is3D)
    
    except IndexError:
        # redo random assemble if the structure 
        # runs out of wild card before includes all fragments
        # e.g. F-F, CF4, CCl4, etc formed
        return RandAssemFrags(fragment_dict, is3D)
        
def CreatePathStateDict(start_mol, end_mol, Gen3D = False):
    start_dict = MakeMolFragDict(PrepMolFrags(start_mol, Gen3D = Gen3D))
    end_dict = MakeMolFragDict(PrepMolFrags(end_mol, Gen3D = Gen3D))

    st_set = set(start_dict.keys())
    end_set = set(end_dict.keys())
    common = st_set.intersection(end_set)
    remove = st_set.difference(end_set)
    add = end_set.difference(st_set)

    common_dict = {}
    remove_dict = {}
    add_dict = {}

    for smi in common:
        mol = start_dict[smi]['mol']
        st_n, end_n = start_dict[smi]['n'], end_dict[smi]['n']
        n_common = min(st_n, end_n)
        common_dict[smi] = {'n': n_common, 'mol': mol}
        n_diff = st_n - end_n
        if n_diff > 0:
            remove_dict[smi] = {'n': n_diff, 'mol': mol}
        elif n_diff < 0:
            add_dict[smi] = {'n': -n_diff, 'mol': mol}

    remove_dict = {
        **remove_dict, 
        **{smi: deepcopy(start_dict[smi]) for smi in remove}
    }

    add_dict = {
        **add_dict, 
        **{smi: deepcopy(end_dict[smi]) for smi in add}
    }
    
    init_state_dict = {
        **deepcopy(start_dict),
        **{smi: {'n': 0, 'mol': end_dict[smi]['mol']} for smi in add}
    }
    
    path_state_dict = {
        'start': start_dict,
        'end': end_dict,
        'common': common_dict,
        'remove': remove_dict,
        'add': add_dict,
        'cur_state': init_state_dict
    }
    
    return path_state_dict

def ReportPathDict(path_state_dict):
    svgs_dict = {}
    for key, frag_dict in path_state_dict.items():
        if frag_dict:
#             print(key)
            svgs_dict[key] = DrawMolFragDict(frag_dict, molsPerRow = 4, useSVG=True)
        else:
#             print(key, 'is EMPTY')
            svgs_dict[key] = None
    return svgs_dict
            
def GenRandPath(perm_path_state_dict, molsPerState = 100, Gen3D = False):
    path_state_dict = deepcopy(perm_path_state_dict)
    choices = []
    for op in ('add','remove'):
        for smi, sub_dict in path_state_dict[op].items():
            for i in range(sub_dict['n']):
                choices.append((op, smi))

    path = {}
    state_id = 0
    while len(choices) > 0:
        op, smi = random.choice(choices)
        choices.remove((op, smi))
        if op == 'add':
            path_state_dict['cur_state'][smi]['n'] += 1
        else:
            path_state_dict['cur_state'][smi]['n'] -= 1

        path_state_dict[op][smi]['n'] -= 1

        if path_state_dict[op][smi]['n'] == 0:
            path_state_dict[op].pop(smi)

        cur_state_mols = []
        for i in range(molsPerState):
            cur_state_mols.append(RandAssemFrags(path_state_dict['cur_state'], is3D = Gen3D))

        path[state_id] = {
            'state': deepcopy(path_state_dict['cur_state']),
            'mols': cur_state_mols
        }
        state_id += 1
    return path