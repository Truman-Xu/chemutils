import os
from rdkit import Chem
from rdkit.Chem.PropertyMol import PropertyMol
from rdkit.Chem import RDConfig # Allow Contrib packages to be used
from rdkit.Chem.Crippen import MolLogP as LogP # Lipophilicity
from rdkit.Chem.QED import default as QED # Quantitiative Estimate of Drug-likeness
from rdkit.Chem.Descriptors import MolWt # Mol Weight
import sys
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
# add path for rdkit Contrib packages
from sascorer import calculateScore as SAS # Sythetic Accessiblilty Score
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map
from itertools import repeat

class PropFilter:
    '''
    Construct a filtering dict with lambda expression for each mol property
    Take input from a dict with filter names (key) followed by a list of relational expressions (value)
    e.g. {'SAS': ["<3"], 'LogP': [">2", "<3"]}
    
    Available Filters: 
        SAS  - Sythetic Accessibility
        QED  - Quantitative Estimate of Druglikeness
        MolWeight - Molecular Weight
        LogP - Octanol-Water Partition Coefficient
    '''
    # Available filters listed below
    # To add a filter, simply give it a name as the key, and the value should be the name of the function for calculating the property
    # The function has to return a numerical value to be able to work in this filtering framework
    filterNames = {
        'MolWeight': MolWt, 
        'LogP': LogP, 
        'QED': QED, 
        'SAS': SAS
    }
    
    def __init__(self, cmd_args_dict):
        self.filter_expr = self._buildFilterDict(cmd_args_dict)
        
    def _validateNum(self, sub_str):
        # Check if the value behind the relational symbols is numeric
        # Also a safefy measure before calling eval() method
        try:
            return float(sub_str)
        except ValueError as e:
            raise ValueError('Filtering criteria contain invalid number: {}\n'.format(sub_str))

    def _assgnBound(self, parsed):   
        # Assign name of the bound based on the symbol
        # return an entry for the dictionary
        if parsed[0] == "<":
            if parsed[1] == "=":
                return ('upper', self._validateNum(parsed[2:]))
            else:
                return ('upper', self._validateNum(parsed[1:]))

        elif parsed[0] == ">":
            if parsed[1] == "=":
                return ('lower', self._validateNum(parsed[2:]))
            else:
                return ('lower', self._validateNum(parsed[1:]))

        else:
            raise ValueError('Invalid Expression: {}. Expression should start with ">", "<", ">=", or "<="\n'.format(parsed))

    def _validateBoundDict(self, bound_dict, args):
        # Check if the bound makes sense
        try:
            return bound_dict['lower'] < bound_dict['upper']
        # error will be raised if two upper or low bounds were specified 
        except KeyError as e:
            raise ValueError('Ambiguous filtering criteria: {}\n'.format(args))

    def _checkAndBuildFilter(self, args):
        # Create a dictionay
        bound_dict = dict(self._assgnBound(parsed) for parsed in args)
        
        if len(args) == 1:
            expression = 'x'+args[0]
        # If two relational expressions specified, check if they make sense
        elif len(args) == 2 and self._validateBoundDict(bound_dict, args):
            expression = ' and '.join('x'+ a for a in args)
        else:
            raise ValueError('\nConflict in lower and upper bound! \nInvalid range: {}'.format(args))
        
        # Create lambda function based on the expression
        return expression

    def _buildFilterDict(self, cmd_args_dict):
        # Create a dict with { filter name : lamdba function to check expression }
        # for filtering each mol property
        expressions = []
        for key, args in cmd_args_dict.items():
            if args and (key in self.filterNames.keys()):
                expr = self._checkAndBuildFilter(args)
                expressions.append((key, expr))
        return dict(expressions)
    
    @staticmethod
    def calcMolProps(mol):
        '''
        Calculate all available filter properties for a molecule
        Does not require instantiation of the class
        Return a dict with all available mol filter properties
        '''
        return dict((key, round(ldFunc(mol),3)) for key, ldFunc in PropFilter.filterNames.items())
    
    def passFilter(self, mol):
        '''
        Compute the properties of the mol and compare to the filter expression
        Return True if all passed and the properties will be assign to the mol object
        Return False if any of the filtering criteria is not met
        '''
        for key in self.filter_expr:
            prop_val = self.filterNames[key](mol)
            ## Multithreading does not work on lambda functions directly
            ## due to the inability to pickle them
            ## eval() method is used as work around
            if not eval("(lambda x: {})(prop_val)".format(self.filter_expr[key])):
                return False
            mol.SetDoubleProp(key, prop_val)
        return True
    
    def filterMols(self, mol_list, show_tqdm = False):
        '''
        Wrapper function for filtering a list of mol based on the built filter dict
        Return a list of mols that pass the filters
        '''
        if show_tqdm:
            return list(filter(self.passFilter, tqdm(mol_list)))
        else:
            return list(filter(self.passFilter, mol_list))

def _sdfFilterWrapper(args):
    propfilter, sdfile_path = args
    mols = [PropertyMol(m) for m in Chem.SDMolSupplier(sdfile_path, removeHs = False) if m]
    return propfilter.filterMols(mols, False)

def sdfFilter(propfilter, sdfile_list, outpath, onefile: bool):
    '''
    Filter molecules in the sdf files based on the criterial specified in propFilter class
    Write the filtered mols in individual sdf files if onefile==False
    Combine the output into a single sdf file if onefile==True
    args:
        -propfilter: the propFilter class built based on the filtering criteria
        -sdfile_list: the list of path to the sdf containing single or multiple mols
        -outpath: the path to the directory to write the filtered mols
        -onefile: if True, combine the all output into one sdf file with the name "filtered.sdf"
    '''
    out_mols = []
    
    if len(sdfile_list) < 30:
        for i, file in enumerate(sdfile_list):
            print('Processing {} of {} sdf files'.format(i+1, len(sdfile_list)))
            cur_mols = [PropertyMol(m) for m in Chem.SDMolSupplier(file, removeHs = False) if m]
            out_mols.append(propfilter.filterMols(cur_mols, show_tqdm = True)) 
                
    else:
        args = list(zip(repeat(propfilter), sdfile_list))
        out_mols = process_map(
            _sdfFilterWrapper,
            args,
            chunksize = 10
        )

        # for i, file in enumerate(tqdm(sdfile_list)):
        #     cur_mols = [m for m in Chem.SDMolSupplier(file, removeHs = False) if m]
        #     out_mols.append(propfilter.filterMols(cur_mols, show_tqdm = False)) 

    if onefile:
        out_file = os.path.join(outpath, 'filtered.sdf')
        print('Writing the combined sdf files')
        with open(out_file, 'w') as f:
            w=Chem.SDWriter(f)
            for ml, sdfile in zip(out_mols, sdfile_list):
                if len(ml) >0:
                    for m in ml:
                        w.write(m)
                else:
                    print(sdfile, "has no qualified mol from the filtering criteria")
            w.close()

    else:
        print('Writing {} sdf files'.format(len(out_mols)))
        for ml, sdfile_path in zip(out_mols, sdfile_list):
            sdfile_name = os.path.split(sdfile_path)[-1]
            if len(ml) >0:
                out_file = os.path.join(outpath, 'filtered_'+sdfile_name)
                with open(out_file, 'w') as f:
                    w=Chem.SDWriter(f)
                    for m in ml:
                        w.write(m)
                    w.close()
            else:
                print(sdfile_name, "has no qualified mol from the filtering criteria")
                
                
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(prog="filterProps",
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                description=
'''
  Physicochemical Filters for Molecules
    Filter molecules in the sdf files based on the criterial specified
    Multiple sdf files are allowed as input
    Write the filtered mols in individual sdf files with prefix "filtered_" on their respective file names.
    Or combine the output into a single sdf file if --onefile is specified

  Available Filters: 
    SAS  - Sythetic Accessibility
    QED  - Quantitative Estimate of Druglikeness
    MW   - Molecular Weight
    LogP - Octanol-Water Partition Coefficient
    
    Filter criteria should start with relational signs, such as ">", "<", ">=", or "<=",
    followed by numerical values, and surrounded by quotation marks
    To specify both upper and lower bound, use two expressions separated by space.
    e.g, for Molecular Weight x: 300 < x < 500, use --MW "<500" ">300"

  Usage Example:
    1.  Process all sdf files in the current directory and write the filtered mols to ./filtered_mols
        filtered by SAS<6, MolWeight<=300, and QED>0.7
        
        python filterProp.py -i ./*.sdf -o ./filtered_mols --SAS "<6" --MW "<=300" --QED ">0.7"
            
    2.  Filter out ligands.sdf and write to current directory as filtered_ligands.sdf (default)
        filtered by MolWeight x: 300 < x < 500, LogP x: 2 <= x <= 3
        
        python filterProp.py -i ligands.sdf --MW ">300" "<500" --LogP "<=3" ">=2"
'''
                                )
    # define extend action for argparse for python version lower than 3.8
    if sys.version_info.major == 3 and sys.version_info.minor < 8:
        class ExtendAction(argparse.Action):
            def __call__(self, parser, namespace, values, option_string=None):
                items = getattr(namespace, self.dest) or []
                items.extend(values)
                setattr(namespace, self.dest, items)
        parser.register('action', 'extend', ExtendAction)

    parser.add_argument('-i','--infile', action="extend", nargs="+", type=str,
                        help='ligand .sdf file path')
    parser.add_argument('-o','--outpath', type=str, default="./",
                        help='directory path for file outputs. Default to current directory')
    parser.add_argument('--onefile', action='store_true', default=False,
                        help='combine output as a single sdf file')
    parser.add_argument('--SAS', '--sas', action="extend", nargs="+", 
                        help='Sythetic Accessibility')
    parser.add_argument('--QED', '--qed', action="extend", nargs="+", type=str,
                        help='Quantitative Estimate of Druglikeness')
    parser.add_argument('--MolWeight', '--MW', '--mw', action="extend", nargs="+", type=str,
                        help='Molecular Weight')
    parser.add_argument('--LogP', '--logP', '--logp', action="extend", nargs="+", type=str,
                        help='Octanol-Water Partition Coefficient')
    
    args = parser.parse_args()
    pf = PropFilter(args.__dict__)
    print('Filters set as below:')
    print(pf.filter_expr)
    
    if not os.path.exists(args.outpath):
        os.makedirs(args.outpath)
        print('Directory Made:',args.outpath)
        
    sdfFilter(pf, args.infile, args.outpath, args.onefile)
    print('Filtering finished!')
    
    