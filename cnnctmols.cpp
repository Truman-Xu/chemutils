#include <iostream>
#include <fstream>
#include <RDGeneral/Invariant.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MolAlign/AlignMolecules.h>
#include <GraphMol/new_canon.h>
#include <GraphMol/FileParsers/MolWriters.h>

using namespace RDKit;

typedef std::map<int, std::pair<int, int>> ORDER_DICT;

void print_map(std::map<int, std::pair<int, int>> m){
    for (const auto& [key, value] : m) {
        std::cout << key << " : " << value.first << ", " << value.second << std::endl;
    }
}

ORDER_DICT findHIds(std::shared_ptr<RDKit::ROMol> mol_ptr, bool unique = true){
    std::vector<unsigned int> rank;
    bool breakTies = !unique;
    Canon::rankMolAtoms(*mol_ptr, rank, breakTies);
    ORDER_DICT orders;
    for (std::size_t j=0; j < rank.size(); j++){
        const Atom *atom = mol_ptr->getAtomWithIdx( j );
        
        if (atom->getSymbol() == "H"){
            // Find the neighbor atom's id
            // There is only one neighbor exists for Hydrogen
            ROMol::ADJ_ITER_PAIR id_pair = mol_ptr->getAtomNeighbors(atom);
            // id of Hydrogen and its neighbor atom matching the sdf idx (H_id, Nbr_id)
            std::pair<int, int> atom_ids{j, *id_pair.first};
            // rank of the hydrogen and the atom id pair
            std::pair<int, std::pair<int, int>> entry{rank[j], atom_ids};
            // Insert in map, only unique rank of hydrogen is allowed as key 
            // Will produce geometrically equivalent hydrogens if unique is false 
            orders.insert(entry);
        }
    }
    return orders;
}

std::shared_ptr<ROMol> connectMols(ROMol &lig, ROMol &frag, 
    int id_ligH, int id_ligNbr, int id_fragH, int id_fragNbr)
    {
        // std::pair<int, int> H_map={id_fragH, id_ligNbr};
        // std::pair<int, int> Nbr_map={id_fragNbr, id_ligH};

        std::shared_ptr<RWMol> rw_lig( new RWMol( lig ));
        std::shared_ptr<RWMol> rw_frag( new RWMol( frag ));

        const MatchVectType id_map = {{id_fragH, id_ligNbr}, {id_fragNbr, id_ligH}};
        double _rmsd = MolAlign::alignMol(*rw_frag, *rw_lig, -1, -1, &id_map);
        
        rw_lig->insertMol( *rw_frag );
        rw_lig->addBond( id_ligNbr, lig.getNumAtoms() + id_fragNbr, Bond::SINGLE );
        rw_lig->removeAtom( lig.getNumAtoms() + id_fragH );
        rw_lig->removeAtom( id_ligH );

        return std::shared_ptr<ROMol>( rw_lig );
}

int main( int argc , char **argv ) {
    std::string lig_file = "/home/truman/chemutils/top10ligs.sdf";
    std::string frag_file = "/home/truman/chemutils/all_frags.sdf";
    bool takeOwnership = true;
    bool removeHs = false;
    SDMolSupplier mol_supplier_lig( lig_file , takeOwnership , removeHs );
    std::shared_ptr<ROMol> lig( mol_supplier_lig.next() );
    
    SDMolSupplier mol_supplier_frag( frag_file , takeOwnership , removeHs );
    std::shared_ptr<ROMol> frag( mol_supplier_frag.next() );
    bool unique=true;
    std::shared_ptr<ROMol> mol2_l( MolOps::removeHs( *lig ));
    std::cout << "lig" << ": " << MolToSmiles(*mol2_l) << std::endl;
    std::shared_ptr<RWMol> mol3_l( new RWMol( *lig ));
    ORDER_DICT orders_lig = findHIds(lig, unique);
    print_map(orders_lig);

    std::shared_ptr<ROMol> mol2_f( MolOps::removeHs( *frag ));
    std::cout << "frag" << ": " << MolToSmiles(*mol2_f) << std::endl;
    std::shared_ptr<RWMol> mol3_f( new RWMol( *frag ));
    ORDER_DICT orders_frag = findHIds(frag, unique);
    print_map(orders_frag);

    std::vector<std::shared_ptr<ROMol>> all_combined;
    for (std::size_t i=0, is=orders_lig.size(); i<is; i++){
        for (std::size_t j=0, js=orders_frag.size(); j<js; j++){

            all_combined.push_back( connectMols(*lig, *frag, 
            orders_lig[i].first, orders_lig[i].second, 
            orders_frag[j].first, orders_frag[j].second) );
        }
    }
    std::string fileout = "combined_out.sdf";
    SDWriter writer(fileout);
    for( std::size_t i = 0 , is = all_combined.size() ; i < is ; ++i ) {
        writer.write( *mols[i] );
    }


    // std::shared_ptr<ROMol> combined( connectMols( *lig, *frag, 
    //     id_ligH, id_ligNbr, id_fragH, id_fragNbr ));

    // std::cout<< lig->getNumAtoms() << std::endl;
    // std::cout<< combined->getNumAtoms() << std::endl;
    
    // std::ofstream ofs( "combined_out.sdf" );
    // ofs << MolToMolBlock( *combined );

    return 0;
}