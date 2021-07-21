#include "functionalize.hpp"
#include <boost/program_options.hpp>

using namespace RDKit;
using namespace chemutils;
namespace po = boost::program_options;

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