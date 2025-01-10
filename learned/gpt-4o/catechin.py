"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    Catechins are flavan-3-ol compounds with specific hydroxyflavan substitution patterns.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """

    # Convert SMILES string to rdkit Molecule object
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for flavan-3-ol skeleton (simplified)
    flavan_3_ol_pattern = Chem.MolFromSmarts("Oc1ccc2c(c1)C[C@@H]3O[C@@H]3c2")
    
    if not mol.HasSubstructMatch(flavan_3_ol_pattern):
        return False, "No flavan-3-ol skeleton match found"

    # Further checks for the presence of typical catechin hydroxylation patterns can be added here
    
    return True, "Matches flavan-3-ol skeleton of catechins"

# Note: While this function checks for the basic flavan-3-ol structure, additional patterns
# and logic would be required to robustly categorize the wide variety of catechin derivatives.