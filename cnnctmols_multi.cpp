#include <iostream>
#include <filesystem>
#include <string>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/FileParsers/MolWriters.h>
#include <thread>
#include <future>
#include <boost/program_options.hpp>

using namespace RDKit;
namespace po = boost::program_options;

typedef std::map<int, std::pair<int, int>> ORDER_DICT;
typedef std::vector<std::shared_ptr<ROMol>> MOL_STDSPTR_VECT;
typedef std::vector<MOL_STDSPTR_VECT> COMBINED_VECT_2D;
typedef std::vector<COMBINED_VECT_2D> COMBINED_VECT_3D;
typedef std::vector<std::shared_future<COMBINED_VECT_2D>> FUTURE_VECT;
typedef std::vector<std::shared_future<void>> FUTURE_WRAPPER_VECT;

ORDER_DICT findHIds(std::shared_ptr<RDKit::ROMol> mol, bool unique = true){
    /*
    Find the hydrogens from the molecule
    If unique flag is set only symmetrically unique hydrogens will be returned
    Return a map with the format {atom order : {H index, Neighbor Atom index}}
    */
    std::vector<unsigned int> rank;
    // find symmetrically equivalent hydrogens by canonical order
    // if unique, no equivalent order will be given (unique rank for equivalent Hs)
    bool breakTies = !unique;
    Canon::rankMolAtoms(*mol, rank, breakTies);
    ORDER_DICT orders;
    for (std::size_t j=0; j < rank.size(); j++){
        const Atom *atom = mol->getAtomWithIdx( j );
        
        if (atom->getSymbol() == "H"){
            // Find the neighbor atom's id
            // There is only one neighbor exists for Hydrogen
            ROMol::ADJ_ITER_PAIR id_pair = mol->getAtomNeighbors(atom);
            // Insert rank of the hydrogen and the atom id pair in map,
            // id of Hydrogen and its neighbor atom matching the sdf idx (H_id, Nbr_id).
            // Only unique rank of hydrogen is allowed as key 
            // Will produce geometrically equivalent hydrogens if unique is false 
            orders.insert({rank[j], {j, *id_pair.first}});
        }
    }
    return orders;
}

void setTetherIdx(std::shared_ptr<ROMol> mol){
    // Set the tethering atom index on the reference lig
    // Mol props set here will be preserved when new is called for constructing 
    // the new ligands
    std::shared_ptr<RDKit::ROMol> mol2( MolOps::removeHs( *mol ));
    // The atom index number to be tethered will be all atoms except hydrogen
    int n_atoms = mol2->getNumAtoms();
    // Construct a string of idx to be set on the mol prop
    std::string tether_idx="1";
    for (int i = 2; i <= n_atoms; i++){
        tether_idx += ", ";
        tether_idx += std::to_string(i);
    }
    // set prop for both versions of rdock
    mol->setProp("rxdock.tethered_atoms", tether_idx);
    mol->setProp("rdock.tethered_atoms", tether_idx);
}

std::shared_ptr<ROMol> connectMols(std::shared_ptr<ROMol> lig, std::shared_ptr<ROMol> frag, 
    int id_ligH, int id_ligNbr, int id_fragH, int id_fragNbr)
    /*
    Connect two mols whose hydrogens are aligned into position
    Return a combined mol
    Got the idea from the python script linked below
    https://github.com/molecularsets/moses/blob/master/moses/baselines/combinatorial.py
    */
    {
        // make editable copied of the molecules
        std::shared_ptr<RWMol> rw_lig( new RWMol( *lig ));
        std::shared_ptr<RWMol> rw_frag( new RWMol( *frag ));
        // map the atom to be aligned from the two mols
        const MatchVectType id_map = {{id_fragH, id_ligNbr}, {id_fragNbr, id_ligH}};
        // mol1->setProp( "rxdock.tethered_atoms", "cyclobutane" );
        // frag first and lig
        double _rmsd = MolAlign::alignMol(*rw_frag, *rw_lig, -1, -1, &id_map);
        // combine the two mols
        rw_lig->insertMol( *rw_frag );
        rw_lig->addBond( id_ligNbr, lig->getNumAtoms() + id_fragNbr, Bond::SINGLE );
        // remove the hydrogens from the ligand and the fragment
        rw_lig->removeAtom( id_fragH + lig->getNumAtoms() );
        rw_lig->removeAtom( id_ligH );
        // sanitize the combined mol 
        try{
            MolOps::sanitizeMol(*rw_lig);
        } catch (MolSanitizeException &e){
            std::cerr << "Fail to sanitize: " << MolToSmiles(*rw_lig) << std::endl;
            std::cerr << e.what() << std::endl;
        }
        
        return std::shared_ptr<ROMol>( rw_lig );
}

MOL_STDSPTR_VECT functionalizeH(std::shared_ptr<ROMol> lig,  
    std::shared_ptr<ROMol> frag, bool unique = true, bool tether=false)
    /*
    Functionalize the ligand with a fragment on all unique positions of hydrogen atoms
    If unqiue flag is set to False, all hydrogen positions on the ligand will be used
    for forming the bonds. 
    Return a vector of ptrs of functionalized ligand
    */
    // For a single lig-frag pair generation
    {
        if (tether == true){
            setTetherIdx(lig);
        }
        ORDER_DICT orders_lig = findHIds(lig, unique);
        // only unique hydrogens position will be used on the fragament
        ORDER_DICT orders_frag = findHIds(frag, true);
        
        MOL_STDSPTR_VECT cur_combined;
        for (ORDER_DICT::iterator it=orders_lig.begin(); it!=orders_lig.end(); ++it){
            for (ORDER_DICT::iterator jt=orders_frag.begin(); jt!=orders_frag.end(); ++jt){
                cur_combined.push_back( connectMols(lig, frag, 
                it->second.first, it->second.second, 
                jt->second.first, jt->second.second) );
            }
        }
        return cur_combined;
}

COMBINED_VECT_2D functionalize_lig(std::shared_ptr<ROMol> lig,  
    MOL_STDSPTR_VECT  frags, bool unique = true, bool tether=false)
    /*
    Enumerated functionalized molecules for a single ligand ptrs 
    and a vector of fragment ptrs
    Return a 2D vector of shared_ptrs with format 
    {lig1-frag1, lig1-frag2, ..., lig1-fragN}
    where each lig-frag pair is a vector of ptrs to functionalized ligands
    */
    {
        if (tether == true){
            setTetherIdx(lig);
        }
        // ligs H ids are determined here to avoid repition
        ORDER_DICT orders_lig = findHIds(lig, unique);
        COMBINED_VECT_2D lig_combined;
        for (std::shared_ptr<ROMol> frag: frags){
            ORDER_DICT orders_frag = findHIds(frag, true);
            MOL_STDSPTR_VECT cur_combined;
            for (ORDER_DICT::iterator it=orders_lig.begin(); it!=orders_lig.end(); ++it){
                for (ORDER_DICT::iterator jt=orders_frag.begin(); jt!=orders_frag.end(); ++jt){
                    cur_combined.push_back( connectMols(lig, frag, 
                    it->second.first, it->second.second, 
                    jt->second.first, jt->second.second) );
                }
            }
            lig_combined.push_back(cur_combined);
        }
        
    return lig_combined;
}

COMBINED_VECT_3D combinatorialFunc(MOL_STDSPTR_VECT ligs, MOL_STDSPTR_VECT frags,
     bool unique = true, bool tether=false){
    /*
    Enumerated all the functionalized molecules 
    from a vector of ligand ptrs and a vector of fragment ptrs
    Return a vector of vectors of shared_ptrs with format 
    {{lig1-frag1, lig1-frag2, ..., lig1-fragN},
     {lig2-frag1, lig2-frag2, ..., lig2-fragN},
      ...
     {ligN-frag1, ligN-frag2, ..., ligN-fragN}}
    */
    int is = ligs.size();
    COMBINED_VECT_3D all_combined;
    FUTURE_VECT fs(is);

    for (int i=0; i<is; i++){
        // multithreading: each thread takes on ligand and functionalize with all fragments
        fs[i] = std::async(functionalize_lig, ligs[i], frags, unique, tether);
    }

    for (int i=0; i<is; i++){
        // collect results
        // threads are async, but the results are placed acording to the ligands ids
        all_combined.push_back(fs[i].get());
    }

    return all_combined;
}

void _sdfCombinatorialWrapper(int lig_idx, std::shared_ptr<ROMol> lig, MOL_STDSPTR_VECT frags, 
    std::string out_path, bool unique, bool tether)
    {
    /*
    Wrapper function to for better multithreading on sdfCombinatorialFunc()
    */
    if (tether == true){
        setTetherIdx(lig);
    }
    ORDER_DICT orders_lig = findHIds(lig, unique);
    int frag_idx = 0;

    for (std::shared_ptr<ROMol> frag: frags){
        // only unique hydrogens position will be used on the fragament
        
        ORDER_DICT orders_frag = findHIds(frag, true);
        std::string fname = "combined_"+std::to_string(lig_idx)+"_"+std::to_string(frag_idx)+".sdf";
        std::string fileout = out_path+"/"+fname;
        SDWriter writer(fileout);

        for (ORDER_DICT::iterator it=orders_lig.begin(); it!=orders_lig.end(); ++it){
            for (ORDER_DICT::iterator jt=orders_frag.begin(); jt!=orders_frag.end(); ++jt){
                writer.write( *connectMols(lig, frag, 
                    // int id_ligH, int id_ligNbr, 
                    // int id_fragH, int id_fragNbr
                    it->second.first, it->second.second, 
                    jt->second.first, jt->second.second)
                );
            }
        }
        frag_idx++;
    }
}

void sdfCombinatorialFunc(std::string lig_file, std::string frag_file, std::string out_path, 
    bool unique=true, bool tether=false){
    /*
    Enumerated all combinatorials of ligands and fragments for various H position for functionalization  
    Takes input from .sdf ligand file and .sdf fragment file 
    Files can contain multiple ligands or fragments to be connected
    Save the fused molecules in sdf files with presevation of ligands 3D coordinates
    */
    if (std::filesystem::exists(out_path)){
        bool takeOwnership = true;
        bool removeHs = false;

        // read ligand mols from file
        MOL_STDSPTR_VECT ligs;
        SDMolSupplier lig_supplier( lig_file , takeOwnership , removeHs );
        while( !lig_supplier.atEnd() ) {
            std::shared_ptr<RDKit::ROMol> lig( lig_supplier.next() );
            if( lig ) {
                ligs.push_back( lig );
            }
        }

        // read fragment mols from file
        MOL_STDSPTR_VECT frags;
        SDMolSupplier frag_supplier( frag_file , takeOwnership , removeHs );
        while( !frag_supplier.atEnd() ) {
            std::shared_ptr<RDKit::ROMol> frag( frag_supplier.next() );
            if( frag ) {
                frags.push_back( frag );
            }
        }
        // execute with multithreaded job
        int is = ligs.size();
        FUTURE_WRAPPER_VECT fs(is);
        for (int i=0; i<is; i++){
            fs[i] = std::async(_sdfCombinatorialWrapper, i, ligs[i], frags, out_path, unique, tether);
        }
        // wait for jobs to finish
        for (int i=0; i<is; i++){
            fs[i].wait();
        }
        std::cout << ligs.size()*frags.size() << " Ligand-Fragament pairs generated." << std::endl;
    } else {
        std::cerr << "[Aborted] File path not exists!" << std::endl;
    }
}

int main( int argc , char **argv ) {
    
    std::string lig_file, frag_file, out_path;
    po::options_description desc("Ligand-Fragment Functionalization");
    bool allH = false, tether = false, verbose=false;
    desc.add_options()
        ("help,h", "produce help message")
        // Option 'buffer-size' and 'b' are equivalent.
        ("lig,l", po::value<std::string>(& lig_file),
            "Ligand sdf to be functionalized.")
        // Option 'priority' and 'p' are equivalent.
        ("frag,f", po::value<std::string>(& frag_file),
            "Fragment sdf to be attached to ligand")
        // Option 'timeout' and 't' are equivalent.
        ("out,o", po::value<std::string>(& out_path)->default_value("./out/"),
            "Directory path for output files. Default to ./out")
        // Option 'allH'
        ("allH", po::bool_switch(& allH),
            "Attach to all Hydrogen position instead of symmetrically unique ones. "
            "Default to false (without this flag)")
        ("tether", po::bool_switch(& tether),
            "Set the tether docking atom ids for the functionalized ligand.\n"
            "If the flag is present, the atoms from the ligand will be fixed during docking, "
            "and only the fragment parts are allowed to move freely.\n"
            "Default to no tethering (without this flag)")
        // Option "verbose". The log messeges need more works
        ("verbose", po::bool_switch(& verbose),
            "Generate more logging messeges"
            "Default to false")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    bool unique= !(allH);
    if (!std::filesystem::exists(out_path)){
        std::filesystem::create_directory(out_path);
        std::cout << "Directory created: "+out_path << std::endl;
    }
    if (verbose) {
        std::cout << "Ligand File: " << lig_file << std::endl;
        std::cout << "Fragment File: " << frag_file << std::endl;
        std::cout << "Out Path: " << out_path << std::endl;
        std::cout << "Attach on All H: " << allH << std::endl;
        std::cout << "Prepare for Tethered Docking: " << tether << std::endl;
    }

    sdfCombinatorialFunc(lig_file, frag_file, out_path, unique, tether);
    return 0;
}