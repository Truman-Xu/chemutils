#ifndef FUNCTIONALIZE_H
#define FUNCTIONALIZE_H

#include <iostream>
#include <string>
#include <filesystem>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <thread>
#include <future>

namespace chemutils
{
typedef std::map<int, std::pair<int, int>> ORDER_DICT;
typedef std::vector<std::shared_ptr<RDKit::ROMol>> MOL_STDSPTR_VECT;
typedef std::vector<MOL_STDSPTR_VECT> COMBINED_VECT_2D;
typedef std::vector<COMBINED_VECT_2D> COMBINED_VECT_3D;
typedef std::vector<std::shared_future<COMBINED_VECT_2D>> FUTURE_VECT;
typedef std::vector<std::shared_future<void>> FUTURE_WRAPPER_VECT;

ORDER_DICT findHIds(std::shared_ptr<RDKit::ROMol> mol, bool unique = true);
/*
    Find the hydrogens from the molecule
    If unique flag is set only symmetrically unique hydrogens will be returned
    Return a map with the format {atom order : {H index, Neighbor Atom index}}
*/

void setTetherIdx(std::shared_ptr<RDKit::ROMol> mol);
/*
    Set the tethering atom index on the reference lig
    Mol props set here will be preserved when new is called for constructing 
    the new ligands
*/

std::shared_ptr<RDKit::ROMol> connectMols(std::shared_ptr<RDKit::ROMol> lig, std::shared_ptr<RDKit::ROMol> frag, 
    int id_ligH, int id_ligNbr, int id_fragH, int id_fragNbr);
/*
    Connect two mols whose hydrogens are aligned into position
    Return a combined mol
    Got the idea from the python script linked below
    https://github.com/molecularsets/moses/blob/master/moses/baselines/combinatorial.py
*/

MOL_STDSPTR_VECT functionalizeH(std::shared_ptr<RDKit::ROMol> lig,  
    std::shared_ptr<RDKit::ROMol> frag, bool unique = true, bool tether=false);
/*
    Functionalize the ligand with a fragment on all unique positions of hydrogen atoms
    If unqiue flag is set to False, all hydrogen positions on the ligand will be used
    for forming the bonds. 
    Return a vector of ptrs of functionalized ligand
    Intended for single lig-frag pair generation
*/

COMBINED_VECT_2D functionalize_lig(std::shared_ptr<RDKit::ROMol> lig,  
    MOL_STDSPTR_VECT  frags, bool unique = true, bool tether=false);
/*
    Enumerated functionalized molecules for a single ligand ptrs 
    and a vector of fragment ptrs
    Return a 2D vector of shared_ptrs with format 
    {lig1-frag1, lig1-frag2, ..., lig1-fragN}
    where each lig-frag pair is a vector of ptrs to functionalized ligands
*/

COMBINED_VECT_3D combinatorialFunc(MOL_STDSPTR_VECT ligs, MOL_STDSPTR_VECT frags,
     bool unique = true, bool tether=false);
/*
    Enumerated all the functionalized molecules 
    from a vector of ligand ptrs and a vector of fragment ptrs
    Return a vector of vectors of shared_ptrs with format 
    {{lig1-frag1, lig1-frag2, ..., lig1-fragN},
     {lig2-frag1, lig2-frag2, ..., lig2-fragN},
      ...
     {ligN-frag1, ligN-frag2, ..., ligN-fragN}}
*/ 

void _sdfCombinatorialWrapper(int lig_idx, std::shared_ptr<RDKit::ROMol> lig, MOL_STDSPTR_VECT frags, 
    std::string out_path, bool unique, bool tether);
/*
    Wrapper function to for better multithreading on sdfCombinatorialFunc()
*/

void sdfCombinatorialFunc(std::string lig_file, std::string frag_file, std::string out_path, 
    bool unique=true, bool tether=false);
    /*
    Enumerated all combinatorials of ligands and fragments for various H position for functionalization  
    Takes input from .sdf ligand file and .sdf fragment file 
    Files can contain multiple ligands or fragments to be connected
    Save the fused molecules in sdf files with presevation of ligands 3D coordinates
    */

} // namespace chemutils

#endif