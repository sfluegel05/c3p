"""
Classifies: CHEBI:26214 porphyrins
"""
from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    Porphyrins have a fundamental skeleton of four pyrrole nuclei connected by four methine groups forming a macrocyclic structure.

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
    
    # Enhanced SMARTS pattern for pyrrole-like units with possible substitutions
    pyrrole_pattern = Chem.MolFromSmarts("c1[c,nH]c[nH]c1")  # More flexible in matching pyrrole variations
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)

    # Ensure at least four pyrrole-like units
    if len(pyrrole_matches) < 4:
        return False, f"Insufficient pyrrole-like units found: {len(pyrrole_matches)} (requires at least 4)"
    
    # SMARTS pattern for methine bridges: a single carbon possibly bonded to other heteroatoms in different porphyrin forms
    methine_bridge_pattern = Chem.MolFromSmarts("C(c1[c,nH]c[nH]c1)C")
    if not mol.HasSubstructMatch(methine_bridge_pattern):
        return False, "Methine bridge pattern not found"
    
    # Verify if there's a macrocyclic structure consisting of the pyrroles and methine bridges
    # Revised logic for macrocycle detection
    cycle_info = mol.GetRingInfo()
    rings = cycle_info.AtomRings()

    # Check for at least one large cycle containing enough atoms (approximately 18-22 typical for porphyrins)
    macrocyclic_detected = any(len(ring) >= 18 and sum(int(cycle_info.IsAtomInRingOfSize(atom, 5)) for atom in ring) >= 4 for ring in rings)
    
    if not macrocyclic_detected:
        return False, "No macrocyclic structure consistent with porphyrin found"
    
    return True, "Contains four pyrrole-like units connected by methine bridges forming a macrocyclic structure"