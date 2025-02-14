"""
Classifies: CHEBI:24128 furanocoumarin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin is a coumarin core with a fused furan ring.

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

    # Define SMARTS patterns for coumarin and furan
    # Coumarin core, allowing for various substitution patterns.
    coumarin_pattern = Chem.MolFromSmarts("[c]1[c][c]([c])[c](=[O])[o][c]1")
    if coumarin_pattern is None:
        return False, "Invalid coumarin SMARTS"
    
    furan_pattern = Chem.MolFromSmarts("[c]1[o][c][c][c]1")  # Furan ring
    if furan_pattern is None:
        return False, "Invalid furan SMARTS"

    # Check for coumarin core
    if not mol.HasSubstructMatch(coumarin_pattern):
        return False, "No coumarin core found"
    
    # Check for furan ring
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found"

    # Check for the fusion. Check for the coumarin and furan with a bond in between.
    #Simplified fusion pattern where there is a bond between carbon of coumarin and carbon of furan
    fused_pattern = Chem.MolFromSmarts("([c]1[c][c]([c])[c](=[O])[o][c]1)-[c]2[o][c][c][c]2")
    if fused_pattern is None:
        return False, "Invalid fused SMARTS"
    
    if not mol.HasSubstructMatch(fused_pattern):
        return False, "Furan ring is not fused to the coumarin core"
    
    return True, "Contains a coumarin core with a fused furan ring"