from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.rdchem import EditableMol, BondType
from rdkit.Chem.MolStandardize import rdMolStandardize


def IsInSameRing(atom_id1, atom_id2, atom_ring_info):
    for t in atom_ring_info:
        if (atom_id1 in t) and (atom_id2 in t):
            return True
    return False

def PrepMolFrags(mol, Gen3D = False, UseWildCard = True, BreakHetSingle = True):
    # Remove Stereochemistry
    Chem.RemoveStereochemistry(mol)
    # Create hydrogen atom as substitution
    subs_atom = Chem.Atom(1) # Hydrogen atom
    mH = AllChem.AddHs(mol)
    rmol = EditableMol(mH)
    atom_ring_info = mol.GetRingInfo().AtomRings()

    if not Gen3D and UseWildCard:
        # Create wild card [*] atom for substitution
        subs_atom = Chem.Atom(0) # Wild card atom
        # Convert all hydrogens to wild card atoms
        for a in mH.GetAtoms():
            if a.GetSymbol() == "H":
                H_idx = a.GetIdx()
                assert len(a.GetNeighbors()) == 1
                Nei_idx = a.GetNeighbors()[0].GetIdx()
                rmol.ReplaceAtom(H_idx, subs_atom)

    # Cleave single bonds and replace with wild card atoms
    # Bonds from the original mol with implicit hydrogens
    # are used here to avoid cleaving the bonds with hydrogen
    for b in mol.GetBonds():
        if b.GetBondType() == BondType.SINGLE:
            bg_idx = b.GetBeginAtomIdx()
            bg_atom = b.GetBeginAtom()
            bg_element = bg_atom.GetSymbol()
            end_idx = b.GetEndAtomIdx()
            end_atom = b.GetEndAtom()
            end_element = end_atom.GetSymbol()
            # Do not break bonds in rings
            if IsInSameRing(bg_idx, end_idx, atom_ring_info):
                continue
            # Do not break single bond between two heteroatoms
            if (bg_element != 'C') and (end_element != 'C'):
                continue
            # Or we want to preserve signle bond between any heteroatoms
            # However, substituent attached to heteroatoms in rings will be cleaved
            if not BreakHetSingle and (
                bg_element != 'C' or end_element != 'C'
            ) and (
                not bg_atom.IsInRing() and not end_atom.IsInRing()
            ):
                continue
            rmol.RemoveBond(bg_idx,end_idx)
            subs_idx_bg = rmol.AddAtom(subs_atom) 
            subs_idx_end = rmol.AddAtom(subs_atom) 
            rmol.AddBond(bg_idx,subs_idx_bg, BondType.SINGLE)
            rmol.AddBond(end_idx,subs_idx_end, BondType.SINGLE)

    out_mol_frags = list(Chem.GetMolFrags(rmol.GetMol(), asMols=True))

    if Gen3D:
        for frag in out_mol_frags:
            AllChem.EmbedMolecule(frag, AllChem.srETKDGv3())

    return out_mol_frags 

def MakeMolFragDict(fragment_list):
    enumerator = rdMolStandardize.TautomerEnumerator()
    frag_dict = {}
    for frag in fragment_list:
        frag = enumerator.Canonicalize(frag)
        smi = Chem.MolToSmiles(frag)
        if smi not in frag_dict:
            frag_dict[smi] = {'n':1, 'mol':frag}
        else:
            frag_dict[smi]['n'] += 1
    return frag_dict

def DrawMolFragDict(frag_dict, add_Hs = False, **kwargs):
    mols, labels = [], []
    for key, sub_dict in frag_dict.items():
        stripped_smi = key.strip('*').replace('(*)','')
        label = f"{stripped_smi}\nn={sub_dict['n']}"
        labels.append(label)
        mol =  Chem.MolFromSmiles(key)
        if add_Hs:
            mol = AllChem.AddHs(mol)
        mols.append(mol)
    return Draw.MolsToGridImage(mols, legends = labels, **kwargs)

def find_largest_ring_size(mol):
    sssr = AllChem.GetSymmSSSR(mol)
    if len(sssr) > 0:
        return max([len(s) for s in sssr])
    return 0

def sort_by_ring_sizes(frag_dict):
    mols = []
    for sub_dict in frag_dict.values():
        mol = sub_dict['mol']
        mol.SetIntProp('occurrence', sub_dict['n'])
        mols.append(mol)

    ring_dict = {
        'none': [],
        'single': [],
        'poly': [],
        'super': []
    }

    for mol in mols:
        ri = mol.GetRingInfo()
        if ri.NumRings() == 0:
            mol.SetIntProp('ring_type', 0)
            ring_dict['none'].append(mol)
        elif ri.NumRings() == 1:
            if find_largest_ring_size(mol) > 8:
                ring_dict['super'].append(mol)
                mol.SetIntProp('ring_type', 3)
            else:
                ring_dict['single'].append(mol)
                mol.SetIntProp('ring_type', 1)
        elif ri.NumRings() > 1:
            if find_largest_ring_size(mol) > 8:
                ring_dict['super'].append(mol)
                mol.SetIntProp('ring_type', 3)
            else:
                ring_dict['poly'].append(mol)
                mol.SetIntProp('ring_type', 2)

    return ring_dict

def has_double_bond(atom):
    for b in atom.GetBonds():
        if b.GetBondType() == BondType.DOUBLE:
            return True
    return False

def get_carbon_nei_group(carbon, group, ter):
    for nei in carbon.GetNeighbors():
        nei_idx = nei.GetIdx()
        if nei.GetSymbol() == 'C':
            if nei_idx not in group:
                group.append(nei_idx)
                group, ter = get_carbon_nei_group(nei, group, ter)
        else:
            ter.append(carbon.GetIdx())
    return group, ter

def find_group_carbons(carbons):
    checked = []
    groups = []
    for carbon in carbons:
        carbon_idx = carbon.GetIdx()
        if carbon_idx in checked:
            continue
        group, ter = get_carbon_nei_group(carbon, [carbon_idx], [])
        groups.append((group, ter))
        checked += group
    return groups


def correct_aromatic_ring(new_mol, check_carbon=True):
    all_pnictogens = []
    protonated = []
    carbons = []
    n_chalcogens = 0
    rings = list(Chem.GetSymmSSSR(new_mol))
    if len(rings) > 1:
        raise ValueError(
                f"Only compound with single ring is supported, {len(rings)} rings found "
            )
    if len(rings) == 0:
        return

    ring_atoms = [new_mol.GetAtomWithIdx(atom_idx) for atom_idx in rings[0]]
    for atom in ring_atoms:
        if not atom.GetIsAromatic():
            return
        symbol = atom.GetSymbol()
        atom.SetFormalCharge(0)
        if symbol == 'N' or symbol == 'P':
            all_pnictogens.append(atom)
            atom.SetNoImplicit(True)
            if atom.GetTotalNumHs() > 0 or len(list(atom.GetNeighbors())) == 3:
                protonated.append(atom)
        elif symbol == 'C':
            if not has_double_bond(atom):
                # if the carbon already has double bond (likely to an oxygen), 
                # no possible aromatic bond contribution in the ring
                carbons.append(atom)
        elif symbol == 'S' or symbol == 'O' or symbol == 'Se':
            n_chalcogens += 1
        else:
            raise ValueError(
                "Only rings containing 'C', 'N', 'O', 'P', 'S', 'Se' are supported "
                f"in correct_aromatic_ring, while {symbol} is found."
            )

    n_pnictogens = len(all_pnictogens)
    n_carbons = len(carbons)
    n_double_bonds = (n_carbons + n_pnictogens) // 2
    n_pnictogens_lone_pairs = 0
    if n_pnictogens:
        n_pnictogens_lone_pairs = (n_carbons + n_pnictogens) % 2
        n_protonation = n_pnictogens_lone_pairs - len(protonated)

        if n_protonation > 0:
            for i in range(abs(n_protonation)):
                all_pnictogens[i].SetNumExplicitHs(1)
        elif n_protonation < 0:
            for i in range(abs(n_protonation)):
                protonated[i].SetNumExplicitHs(0)

    if check_carbon and n_carbons > 0:
        carbon_groups = find_group_carbons(carbons)
        for group, ter in carbon_groups:
            if len(group) == 1:
                selected_carbon = new_mol.GetAtomWithIdx(group[0])
                selected_carbon.SetFormalCharge(2)
            elif len(group) % 2:
                if len(ter) > 0:
                    selected_carbon = new_mol.GetAtomWithIdx(ter[0])
                else:
                    selected_carbon = new_mol.GetAtomWithIdx(group[0])
                selected_carbon.SetFormalCharge(2)

    n_total_lone_pairs = n_pnictogens_lone_pairs + n_chalcogens
    n_total_electrons = (n_total_lone_pairs + n_double_bonds) * 2
    return ((n_total_electrons - 2) % 4 == 0)

def remove_aromaticity(mol):
    for atom in mol.GetAtoms():
        atom.SetIsAromatic(False)
    for bond in mol.GetBonds():
        if bond.GetBondType() == BondType.AROMATIC:
            bond.SetBondType(BondType.SINGLE)

def find_ring_atom_neighbor(atom, checked):

    atom_idx = atom.GetIdx()
    if atom.IsInRing():
        return atom_idx
    checked.append(atom_idx)

    neighbors = list(atom.GetNeighbors())
    if len(neighbors) == 0:
        return

    for nei in neighbors:
        if nei.GetIdx() in checked:
            continue
        return find_ring_atom_neighbor(nei, checked)

def find_non_ring_atom_roots(mol):
    non_ring_atoms = {}
    for atom in mol.GetAtoms():
        if not atom.IsInRing():
            cur_group = []
            ring_atom_idx = find_ring_atom_neighbor(atom, cur_group)
            if ring_atom_idx:
                non_ring_atoms[ring_atom_idx] = cur_group
    return non_ring_atoms

def sep_poly_ring_atom_ids(mol):
    roots = find_non_ring_atom_roots(mol)
    rings = []
    for ring in Chem.GetSymmSSSR(mol):
        ring = set(ring)
        for root, substituents in roots.items():
            if root in ring:
                ring = ring.union(substituents)
        rings.append(sorted(ring))
    return rings

def sep_poly_rings(mol, check_carbon=True):
    sub_mols = []
    dirty_mols = []
    for sub_mol_atom_ids in sep_poly_ring_atom_ids(mol):
        smiles = Chem.MolFragmentToSmiles(mol, sub_mol_atom_ids, kekuleSmiles=False)
        sub_mol = Chem.MolFromSmiles(smiles, sanitize=False)
        error_code = Chem.SanitizeMol(sub_mol, catchErrors=True)
        if error_code:
            sub_mol = Chem.MolFromSmiles(smiles, sanitize=False)
            is_aro = correct_aromatic_ring(sub_mol, check_carbon=check_carbon)
            if is_aro is None:
                remove_aromaticity(sub_mol)

            error_code = Chem.SanitizeMol(sub_mol, catchErrors=True)
        if error_code:
            dirty_mols.append(sub_mol)
        else:
            sub_mols.append(sub_mol)
    return sub_mols, dirty_mols

def _has_carbon(mol):
    """Check if the molecule contains carbon atom"""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            return True
    return False

def filter_inorganic_mols(mol_list):
    """Filter out inorganic molecules from a list of molecules. Single atom
    molecules are also removed."""
    filtered = []
    for mol in mol_list:
        if mol.GetNumAtoms() == 1:
            continue
        elif not _has_carbon(mol):
            continue
        filtered.append(mol)
    return filtered
