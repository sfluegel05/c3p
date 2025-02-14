"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is defined as a six-carbon monosaccharide with an aldehyde at position 1 (aldohexose)
    or a ketone group at position 2 (ketohexose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check if there are exactly six carbon atoms
    if c_count != 6:
        return False, f"Contains {c_count} carbon atoms, but a hexose requires exactly 6"

    # Define patterns for an aldehyde at position 1 and a ketone at position 2
    aldehyde_pattern = Chem.MolFromSmarts("C(=O)[CH]")
    ketone_pattern = Chem.MolFromSmarts("C[C](=O)C")

    # Check for an aldehyde at position 1
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Structure matches an aldohexose (aldehyde at position 1)"
    
    # Check for a ketone at position 2
    if mol.HasSubstructMatch(ketone_pattern):
        return True, "Structure matches a ketohexose (ketone at position 2)"

    return False, "Does not contain hexose-defining functional groups (aldehyde at C1 or ketone at C2)"