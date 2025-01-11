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
    
    # Define SMARTS pattern for a furanocoumarin backbone
    # Patterns should capture the various potential fusion types of furo and coumarin
    furanocoumarin_patterns = [
        Chem.MolFromSmarts("c1oc2ccccc2oc1C=O"), # Basic furo[2,3-b]coumarin structure
        Chem.MolFromSmarts("c1oc2cccc(c2cc1C=O)"), # Basic furo[3,2-g]coumarin structure
        Chem.MolFromSmarts("c1oc2ccc(o2)c(c1)C=O")  # Basic psoralen-like structure
    ]
    
    # Ensure all patterns are valid
    if any(pattern is None for pattern in furanocoumarin_patterns):
        return None, None  # If any pattern is None, there's a problem with its definition
    
    # Check for presence of any furanocoumarin patterns
    for pattern in furanocoumarin_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains furanocoumarin fused rings"
            
    return False, "The structure lacks necessary furanocoumarin fused moieties"