"""
Classifies: CHEBI:26214 porphyrins
"""
from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    Porphyrins have a fundamental skeleton of four pyrrole nuclei connected by methine bridges, typically forming a macrocyclic structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Update pyrrole-like pattern allowing more flexibility
    pyrrole_pattern = Chem.MolFromSmarts("[n;R1]1[c][c][c][c]1")
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)

    if len(pyrrole_matches) < 4:
        return False, f"Insufficient pyrrole-like units found: {len(pyrrole_matches)} (requires at least 4)"
    
    # Update methine bridge pattern for exact porphyrin connector: Methine bridges must be sp2 carbon connected to pyrroles
    methine_pattern = Chem.MolFromSmarts("[C]=[C]([n])[C]")  # Capture the connecting carbons more specifically
    methine_matches = mol.GetSubstructMatches(methine_pattern)

    if len(methine_matches) < 4:
        return False, f"Insufficient methine bridge features found: {len(methine_matches)} (requires at least 4)"

    # Check for macrocyclic structure using ring detection methods
    cycle_info = mol.GetRingInfo()
    rings = cycle_info.AtomRings()

    macrocyclic_detected = any(
        len(ring) >= 16 and
        sum(int(cycle_info.IsAtomInRingOfSize(atom, 5)) for atom in ring) >= 4
        for ring in rings
    )

    if not macrocyclic_detected:
        return False, "No macrocyclic structure consistent with porphyrin found"
    
    return True, "Contains four pyrrole-like units connected by methine bridges forming a macrocyclic structure"