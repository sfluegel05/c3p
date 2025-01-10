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

    # Get the ring information to ensure presence of 5 or 6 member rings
    ring_info = mol.GetRingInfo()
    rings = [set(indices) for indices in ring_info.AtomRings() if len(indices) in [5, 6]]
    if not rings:
        return False, "No suitable 5 or 6 membered ring structure found"

    # Define a SMARTS pattern for lactol structure
    # Oxygen inside a ring connected to an alcohol and a secondary carbon
    lactol_pattern = Chem.MolFromSmarts("OC1(O)[C@H1]")

    # Check if the molecule matches the lactol pattern in an appropriate ring
    if not any(mol.HasSubstructMatch(lactol_pattern, atoms=ring) for ring in rings):
        return False, "No lactol pattern found (cyclic hemiacetal)"

    return True, "Contains a lactol structural pattern (cyclic hemiacetal with adjacent hydroxyl) in the molecule"