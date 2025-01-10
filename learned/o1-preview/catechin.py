"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: catechin
"""
from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    A catechin is a hydroxyflavan with a flavan-3-ol skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for flavan-3-ol skeleton without stereochemistry
    # This pattern represents the core flavan-3-ol structure
    flavan_3_ol_smarts = 'c1ccc2c(c1)CC(O)CO2'
    flavan_3_ol_pattern = Chem.MolFromSmarts(flavan_3_ol_smarts)
    if flavan_3_ol_pattern is None:
        return False, "Invalid SMARTS pattern for flavan-3-ol skeleton"
    
    # Check if molecule has the flavan-3-ol substructure
    if not mol.HasSubstructMatch(flavan_3_ol_pattern):
        return False, "Flavan-3-ol skeleton not found"
    
    return True, "Molecule contains flavan-3-ol skeleton characteristic of catechins"

# Example usage:
# smiles_example = 'O[C@H]1Cc2c(O)cc(O)cc2O[C@@H]1c1ccc(O)c(O)c1'  # (+)-catechin
# result = is_catechin(smiles_example)
# print(result)