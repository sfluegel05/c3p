"""
Classifies: CHEBI:71971 neoflavonoid
"""
from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a coumarin-like pattern (primary structure for many neoflavonoids)
    coumarin_pattern = Chem.MolFromSmarts('O=C1OC2=CC=CC=C2C1=O')
    if not mol.HasSubstructMatch(coumarin_pattern):
        return False, "No coumarin-like backbone found"

    # Look for phenyl group (aromatic ring) which is a common substitution on coumarins
    phenyl_pattern = Chem.MolFromSmarts('c1ccccc1')
    if not mol.HasSubstructMatch(phenyl_pattern):
        return False, "No phenyl group substitution found"

    return True, "Contains coumarin-like structure with phenyl substitution"

# Examples: The provided examples can be tested to verify the functionality.