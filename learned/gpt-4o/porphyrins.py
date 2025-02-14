"""
Classifies: CHEBI:26214 porphyrins
"""
from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    Porphyrins have a fundamental skeleton of four pyrrole-like units connected by four methine groups forming a macrocyclic structure.

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

    # Define pyrrole-like pattern more flexibly
    pyrrole_pattern = Chem.MolFromSmarts("[nH]1[cR1][cR1][cR1][cR1]1")  # Allowing for aromatic conjugation
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)

    if len(pyrrole_matches) < 4:
        return False, f"Insufficient pyrrole-like units found: {len(pyrrole_matches)} (requires at least 4)"
  
    # Methine bridge pattern: Carbon bonded to other carbon/nitrogen that forms part of the macrocycle
    methine_bridge_pattern = Chem.MolFromSmarts("C~[C,N]")  # More generic linkage
    methine_matches = mol.GetSubstructMatches(methine_bridge_pattern)

    # Check for minimum number of methine bridges typically seen in porphyrins
    if len(methine_matches) < 4:
        return False, f"Insufficient methine bridge-like features found: {len(methine_matches)} (requires at least 4)"

    # Enhanced cycle detection logic for macrocyclic structure comprised of identified units
    cycle_info = mol.GetRingInfo()
    rings = cycle_info.AtomRings()

    # Check for macrocycle with enough pyrrole-like units (using typical porphyrin ring sizes)
    macrocyclic_detected = any(
        len(ring) >= 18 and 
        sum(int(cycle_info.IsAtomInRingOfSize(atom, 5)) for atom in ring) >= 4
        for ring in rings
    )
    
    if not macrocyclic_detected:
        return False, "No macrocyclic structure consistent with porphyrin found"
    
    return True, "Contains four pyrrole-like units connected by methine bridges forming a macrocyclic structure"