#include <iostream>
// #include <fstream>
#include <filesystem>
// #include <RDGeneral/Invariant.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
// #include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/MolSupplier.h>
// #include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/FileParsers/MolWriters.h>

using namespace RDKit;

typedef std::map<int, std::pair<int, int>> ORDER_DICT;
typedef std::vector<std::shared_ptr<ROMol>> MOL_VEC;

void print_map(std::map<int, std::pair<int, int>> m){
    for (const auto& [key, value] : m) {
        std::cout << key << " : " << value.first << ", " << value.second << std::endl;
    }
}

ORDER_DICT findHIds(std::shared_ptr<RDKit::ROMol> mol, bool unique = true){
    std::vector<unsigned int> rank;
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

std::shared_ptr<ROMol> connectMols(std::shared_ptr<ROMol> lig, std::shared_ptr<ROMol> frag, 
    int id_ligH, int id_ligNbr, int id_fragH, int id_fragNbr)
    
    {
        // make editable copied of the molecules
        std::shared_ptr<RWMol> rw_lig( new RWMol( *lig ));
        std::shared_ptr<RWMol> rw_frag( new RWMol( *frag ));
        // map the atom to be aligned from the two mols
        const MatchVectType id_map = {{id_fragH, id_ligNbr}, {id_fragNbr, id_ligH}};
        // frag first and lig
        double _rmsd = MolAlign::alignMol(*rw_frag, *rw_lig, -1, -1, &id_map);
        
        rw_lig->insertMol( *rw_frag );
        rw_lig->addBond( id_ligNbr, lig->getNumAtoms() + id_fragNbr, Bond::SINGLE );
        rw_lig->removeAtom( id_fragH + lig->getNumAtoms() );
        rw_lig->removeAtom( id_ligH );
        return std::shared_ptr<ROMol>( rw_lig );
}

MOL_VEC functionalizeH(std::shared_ptr<ROMol> lig,  
    std::shared_ptr<ROMol> frag, bool unique = true)
    /*
    Functionalize the ligand with a fragment on all unique positions of hydrogen atoms
    If unqiue flag is set to False, all hydrogen positions on the ligand will be used
    for forming the bonds. 
    Return a vector of ptrs of functionalized ligand
    */
    {
        ORDER_DICT orders_lig = findHIds(lig, unique);
        ORDER_DICT orders_frag = findHIds(frag, true);

        MOL_VEC all_combined;
        for (ORDER_DICT::iterator it=orders_lig.begin(); it!=orders_lig.end(); ++it){
            for (ORDER_DICT::iterator jt=orders_frag.begin(); jt!=orders_frag.end(); ++jt){
                all_combined.push_back( connectMols(lig, frag, 
                it->second.first, it->second.second, 
                jt->second.first, jt->second.second) );
            }
        }
        return all_combined;
}

void sdfEnumFuncH(std::string lig_file, std::string frag_file, std::string out_path, bool unique=true){
    if (std::filesystem::exists(out_path)){
        bool takeOwnership = true;
        bool removeHs = false;
        MOL_VEC ligs;
        SDMolSupplier lig_supplier( lig_file , takeOwnership , removeHs );
        while( !lig_supplier.atEnd() ) {
            std::shared_ptr<RDKit::ROMol> lig( lig_supplier.next() );
            if( lig ) {
                ligs.push_back( lig );
            }
        }

        MOL_VEC frags;
        SDMolSupplier frag_supplier( frag_file , takeOwnership , removeHs );
        while( !frag_supplier.atEnd() ) {
            std::shared_ptr<RDKit::ROMol> frag( frag_supplier.next() );
            if( frag ) {
                frags.push_back( frag );
            }
        }

        for (std::size_t i=0, is=ligs.size(); i<is; i++){
            for (std::size_t j=0, js=frags.size(); j<js; j++){
                MOL_VEC combined = functionalizeH( ligs[i], frags[j], unique );
                std::string fname = "combined_"+std::to_string(i)+"_"+std::to_string(j)+".sdf";
                std::string fileout = out_path+"/"+fname;
                SDWriter writer(fileout);
                for( std::size_t k = 0 , ks = combined.size() ; k < ks ; ++k ) {
                    writer.write( *combined[k] );
                }
            }
        }
        std::cout << ligs.size()*frags.size() << " Ligand-Fragament pairs generated." << std::endl;
    } else {
        std::cout << "[Aborted] File path not exists!" << std::endl;
    }
}

int main( int argc , char **argv ) {
    for (int i = 0; i < argc; ++i)
        std::cout << argv[i] << "\n";

    std::string lig_file = argv[1], frag_file = argv[2], out_path = argv[3];

    bool unique=true;
    std::string isUnique;
    if(argc>3){
        isUnique = argv[4];
        if (isUnique == "true" || isUnique == "True") { 
            unique = true;
        } else if (isUnique == "false" || isUnique == "False") {
            unique = false;
        } else {
            // Report invalid argument
            std::cout << "Invalid input for unique Hydrogen, only true or false is allowed." << std::endl;
            std::cout << "Using default value: true" << std::endl;
        }
    }
    
    if (!std::filesystem::exists(out_path)){
        std::filesystem::create_directory(out_path);
        std::cout << "Directory created: "+out_path << std::endl;
    }
    
    sdfEnumFuncH(lig_file, frag_file, out_path, unique);
    
    return 0;
}