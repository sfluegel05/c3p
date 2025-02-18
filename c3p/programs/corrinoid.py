"""
Classifies: CHEBI:33913 corrinoid
"""
from rdkit import Chem

def is_corrinoid(smiles: str):
    """
    Determines if a molecule is a corrinoid based on its SMILES string.
    A corrinoid typically contains a corrin macrocycle with cobalt.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a corrinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Corrin ring pattern is hard to capture in a singular SMARTS
    # We can check for cobalt and a large macrocycle as proxy
    cobalt_atom = None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'Co':
            cobalt_atom = atom
            break
    if not cobalt_atom:
        return False, "Cobalt atom not found"

    # Attempt to identify a macrocyclic ring structure (indicative of a corrin)
    # Corrin-specific patterns would make use of detailed SMARTS capturing the macrocycle
    # This is a simplified check to capture a large ring presence
    ri = mol.GetRingInfo()
    if not any(len(r) > 10 for r in ri.AtomRings()):
        return False, "No large macrocyclic ring found"

    return True, "Contains cobalt and large macrocyclic ring consistent with corrinoid structure"