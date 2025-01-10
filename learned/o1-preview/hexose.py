"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is any six-carbon monosaccharide which in its linear form contains either
    an aldehyde group at position 1 (aldohexose) or a ketone group at position 2 (ketohexose).

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

    # Count number of carbons
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons != 6:
        return False, f"Number of carbon atoms ({num_carbons}) is not 6"

    # Check for aldehyde group at position 1 (terminal aldehyde)
    aldehyde_pattern = Chem.MolFromSmarts("[#6H1][CX3H1](=O)")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if aldehyde_matches:
        # Confirm aldehyde is at position 1
        return True, "Molecule is an aldohexose (contains terminal aldehyde group at position 1)"

    # Check for ketone group at position 2
    ketose_pattern = Chem.MolFromSmarts("[#6H2][#6](=O)[#6H1]")
    ketose_matches = mol.GetSubstructMatches(ketose_pattern)
    if ketose_matches:
        # Confirm ketone is at position 2
        return True, "Molecule is a ketohexose (contains ketone group at position 2)"

    return False, "Molecule does not have an aldehyde at position 1 or ketone at position 2"