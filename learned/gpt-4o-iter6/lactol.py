"""
Classifies: CHEBI:38131 lactol
"""
from rdkit import Chem

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
    if not ring_info.IsAtomInRingOfSize(0, 5) and not ring_info.IsAtomInRingOfSize(0, 6):
        return False, "No suitable ring structure found - lactols must be 5 or 6 member cyclic"

    # Define a SMARTS pattern for lactol structure
    # This pattern detects a cyclic ether linked to a hydroxyl group
    # within 5 or 6 membered rings which are typical for lactols.
    lactol_pattern = Chem.MolFromSmarts("O[C]1[OH][C,O][c,C]1")

    # Check if the molecule matches the lactol pattern
    if not mol.HasSubstructMatch(lactol_pattern):
        return False, "No lactol pattern found (cyclic hemiacetal)"

    return True, "Contains a lactol structural pattern (cyclic hemiacetal with adjacent hydroxyl) in the molecule"