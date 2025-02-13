"""
Classifies: CHEBI:50523 butenolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is defined as a gamma-lactone that consists of a 2-furanone skeleton
    and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define gamma-lactone (lactone) pattern (O=C1OC(C=C1))
    lactone_pattern = Chem.MolFromSmarts("O=C1OC=C1")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No gamma-lactone (lactone) structure found"
    
    # Define 2-furanone skeleton pattern (O=C1OC=CC1)
    furanone_pattern = Chem.MolFromSmarts("O=C1OC=CC1")
    if not mol.HasSubstructMatch(furanone_pattern):
        return False, "No 2-furanone skeleton found"
    
    # If both patterns are found, it is a butenolide
    return True, "Contains gamma-lactone structure with 2-furanone skeleton"

# Example usage
smiles = "O=C1OCC=C1"  # Example of butenolide (but-2-en-4-olide)
is_butenolide(smiles)