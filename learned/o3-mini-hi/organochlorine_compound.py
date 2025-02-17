"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies: organochlorine compound
Definition: An organochlorine compound is one that contains at least one carbon-chlorine bond.
"""
from rdkit import Chem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule contains at least one carbon-chlorine bond, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Create a SMARTS pattern that matches a C-Cl bond.
    # [#6] matches any atom with atomic number 6 (carbon).
    # -Cl ensures there is a single bond connecting to a chlorine atom.
    ccl_pattern = Chem.MolFromSmarts("[#6]-Cl")
    
    # Check whether the molecule contains at least one carbon-chlorine bond.
    if mol.HasSubstructMatch(ccl_pattern):
        return True, "Contains at least one carbon-chlorine bond"
    else:
        return False, "Does not contain any carbon-chlorine bonds"
    
# Example usage (uncomment below to test):
# test_smiles = "ClCCBr"  # 1-bromo-2-chloroethane example
# print(is_organochlorine_compound(test_smiles))