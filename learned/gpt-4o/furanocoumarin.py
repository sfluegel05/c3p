"""
Classifies: CHEBI:24128 furanocoumarin
"""
from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin is defined as a furochromene that consists of a furan ring 
    fused with a coumarin.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecular object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Expanded SMARTS pattern for a furanocoumarin backbone
    # Covering a broader range of furan and coumarin linkages
    furanocoumarin_patterns = [
        Chem.MolFromSmarts("c1oc2ccc3c(c2oc1)C=CC(=O)c3"), # General fused furo[3,2-g]chromene
        Chem.MolFromSmarts("c1cc2ccoc2c3c1ccco3"),         # Furo[2,3-b]chromene structure
        Chem.MolFromSmarts("c1occ2c(c1)ccc3ccoc23"),       # Psoralen-like and other isomeric structures
        Chem.MolFromSmarts("c1oc2cc3c(cc2o1)c(=O)cco3")    # Additional fusion pattern
    ]
    
    # Ensure all patterns are valid
    if any(pattern is None for pattern in furanocoumarin_patterns):
        return None, None  # If any pattern is None, there's a problem with its definition
    
    # Check for presence of any furanocoumarin patterns
    for pattern in furanocoumarin_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains furanocoumarin fused rings"
            
    return False, "The structure lacks necessary furanocoumarin fused moieties"