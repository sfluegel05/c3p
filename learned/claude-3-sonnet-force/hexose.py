"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is a six-carbon monosaccharide with either an aldehyde group at position 1 (aldohexose)
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
        return False, "Invalid SMILES string"

    # Check for 6 carbon atoms and at least one oxygen
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count != 6 or o_count < 1:
        return False, "Not a six-carbon monosaccharide"

    # Check for aldehyde or ketone group
    aldehyde_pattern = Chem.MolFromSmarts("C=O")
    ketone_pattern = Chem.MolFromSmarts("[CX3]=O")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)

    if len(aldehyde_matches) + len(ketone_matches) == 0:
        return False, "No aldehyde or ketone group found"

    # Check for linear form (optional)
    # ...

    # Check for specific ring structures (optional)
    # ...

    # If all conditions are met, classify as hexose
    if len(aldehyde_matches) > 0:
        return True, "Molecule is an aldohexose"
    else:
        return True, "Molecule is a ketohexose"