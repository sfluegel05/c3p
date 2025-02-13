"""
Classifies: CHEBI:38131 lactol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is a cyclic hemiacetal formed by intramolecular addition of a hydroxy group
    to an aldehydic or ketonic carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ring systems
    ring_info = mol.GetRingInfo()
    if not ring_info:
        return False, "No ring structure found - lactols must be cyclic"

    # Define a Smart pattern for a lactol group
    # This pattern must detect rings with an oxygen atom within the ring forming a cyclic ether
    # and an adjacent hydroxyl group indicative of the lactol structure.
    lactol_pattern = Chem.MolFromSmarts("O1[C@@H](O)C1")

    # Check if the molecule matches the lactol pattern
    if not mol.HasSubstructMatch(lactol_pattern):
        return False, "No lactol pattern found (cyclic hemiacetal)"

    return True, "Contains a lactol structural pattern (cyclic hemiacetal with adjacent hydroxyl) in the molecule"