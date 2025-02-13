"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: Furanocoumarin
"""
from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin consists of a furan ring fused with a coumarin structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Coumarin structure pattern (benzopyrone structure)
    coumarin_pattern = Chem.MolFromSmarts("O=c1oc2ccccc2cc1")

    # Furan heterocyclic ring
    furan_pattern = Chem.MolFromSmarts("c1occc1")
    
    # Patterns for linear and angular furanocoumarin
    psoralen_pattern = Chem.MolFromSmarts("O=c1oc2cc3c(ccc3oc2c4ccccc4)c1")  # Linear type
    angelicin_pattern = Chem.MolFromSmarts("O=c1c2c(ccc3oc4c(c2)ccoc4cc3)o1") # Angular type
    
    # Check for basic coumarin structure
    if not mol.HasSubstructMatch(coumarin_pattern):
        return False, "No coumarin structure found"
    
    # Check for general furan ring
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found"

    # Check for linear or angular furanocoumarin pattern
    if mol.HasSubstructMatch(psoralen_pattern):
        return True, "Contains linear furanocoumarin (psoralen) structure"
    
    if mol.HasSubstructMatch(angelicin_pattern):
        return True, "Contains angular furanocoumarin (angelicin) structure"

    return False, "No recognized furanocoumarin fusion detected"