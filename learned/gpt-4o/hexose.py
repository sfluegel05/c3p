"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is defined as a six-carbon monosaccharide which in its linear form
    contains either an aldehyde group at position 1 (aldohexose) or a ketone group at position 2 (ketohexose).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, f"Expected 6 carbon atoms, found {c_count}"

    # Define aldehyde and ketone patterns
    aldehyde_pattern = Chem.MolFromSmarts("[#6H1](=O)")
    ketone_pattern = Chem.MolFromSmarts("[#6H0,C]([#6H0,C])=O")

    # Check for either aldehyde or ketone group
    is_aldohexose = mol.HasSubstructMatch(aldehyde_pattern)
    is_ketohexose = mol.HasSubstructMatch(ketone_pattern)

    if is_aldohexose:
        return True, "Contains aldohexose group"
    elif is_ketohexose:
        return True, "Contains ketohexose group"
    else:
        return False, "Does not contain aldehyde or ketone group in the right position"