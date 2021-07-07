#include <iostream>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/new_canon.h>

using namespace RDKit;

typedef std::map<int, std::pair<int, int>> ORDER_DICT;

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
    for (int j=0; j<rank.size(); j++){
        const Atom *atom = mol->getAtomWithIdx( j );
        
        if (atom->getSymbol() == "H"){
            // Find the neighbor atom's id
            // There is only one neighbor exists for Hydrogen
            ROMol::ADJ_ITER_PAIR id_pair = mol->getAtomNeighbors(atom);
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

int main( int argc , char **argv ) {
    // Check the number of parameters
    if (argc < 2) {
        // Tell the user how to run the program
        std::cerr << "Usage: " << argv[0] << " ligand sdf file" << std::endl;
        /* "Usage messages" are a conventional way of telling the user
         * how to run a program if they enter the command incorrectly.
         */
        return 1;
    }

    std::string sdf_file = argv[1];
    bool unique;
    std::string isUnique(argv[2]);
    if (isUnique == "true" || isUnique == "True") { 
        unique = true;
    } else if (isUnique == "false" || isUnique == "False") {
        unique = false;
    } else {
        // Report invalid argument
        std::cout << "Invalid input for unique Hydrogen, only true or false is allowed." << std::endl;
        std::cout << "Using default value: true" << std::endl;
    }

    bool takeOwnership = true;
    bool removeHs = false;
    SDMolSupplier mol_supplier( sdf_file , takeOwnership , removeHs );

    std::vector<std::shared_ptr<RDKit::ROMol>> mols;
    while( !mol_supplier.atEnd() ) {
        std::shared_ptr<RDKit::ROMol> mol( mol_supplier.next() );
        if( mol ) {
            mols.push_back( mol );
        }
    }
    
    for( std::size_t i = 0 , is = mols.size() ; i < is ; ++i ) {
        std::shared_ptr<ROMol> mol2( MolOps::removeHs( *mols[i] ));
        std::cout << i << ": " << MolToSmiles(*mol2) << std::endl;
        ORDER_DICT orders = findHIds(mols[i], unique);
        print_map(orders);
        // if (atom->getSymbol() == "H"){
            // const Atom *atom = mol->getAtomWithIdx( j );
            // ptr to the neighbor atom
            // const Atom *nbr = (*mols[i])[*id_pair.first];
            // std::cout << ' ' << rank[j] << ": " << atom->getSymbol() << ", " 
            // << nbr->getSymbol() << " (" << j  << ", " << *id_pair.first  << ")" << std::endl;
        // }
    }
    return 0;
}