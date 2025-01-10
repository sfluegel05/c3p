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
    # Start with capturing typical furan rings fused to coumarin in several potential configurations
    furanocoumarin_patterns = [
        Chem.MolFromSmarts("c1c2occc2c(c2oc(=O)cc)cc1"),  # Pattern for psoralen-like structure
        Chem.MolFromSmarts("c1oc2ccc(o2)c(c1)c(=O)c"),    # Pattern for a basic furochromene
    ]
    
    # Check for presence of any furanocoumarin patterns
    for pattern in furanocoumarin_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains furanocoumarin fused rings"
            
    return False, "The structure lacks necessary furanocoumarin fused moieties"