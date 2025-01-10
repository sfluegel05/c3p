"""
Classifies: CHEBI:24128 furanocoumarin
"""
from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin consists of a furan ring fused with a coumarin, allowing for different fusion patterns.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for different types of furanocoumarin structures
    psoralen_pattern = Chem.MolFromSmarts('c1cc2oc(=O)ccc2oc1')  # Linear furanocoumarins
    angelicin_pattern = Chem.MolFromSmarts('c1cc2oc(=O)c3c1ccco3c2')  # Angular furanocoumarins

    # Check for the presence of any recognized furanocoumarin core structure
    if mol.HasSubstructMatch(psoralen_pattern):
        return True, "Contains furanocoumarin core (linear furanocoumarin)"
    if mol.HasSubstructMatch(angelicin_pattern):
        return True, "Contains furanocoumarin core (angular furanocoumarin)"
    
    return False, "No furanocoumarin core structure found"

# This function now accounts for both major furanocoumarin core types, improving its ability to classify various isomers and derivatives.