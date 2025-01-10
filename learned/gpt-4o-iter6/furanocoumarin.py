"""
Classifies: CHEBI:24128 furanocoumarin
"""
from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin is defined as a furan ring fused with a coumarin structure.

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

    # Look for furan ring pattern (5-membered oxygen-containing ring)
    furan_pattern = Chem.MolFromSmarts("o1cccc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found"
    
    # Look for coumarin structure pattern (benzopyran-2-one)
    coumarin_pattern = Chem.MolFromSmarts("O=c1ccc2c(c1)occ2")
    if not mol.HasSubstructMatch(coumarin_pattern):
        return False, "No coumarin structure found"
    
    # Check for fusion of furan and coumarin
    fused_pattern = Chem.MolFromSmarts("o1c2c(cccc2)oc(=O)c3c1ccc2c3cccc2")
    if not mol.HasSubstructMatch(fused_pattern):
        return False, "Furan and coumarin are not correctly fused"

    return True, "Contains a furan ring fused with a coumarin structure"