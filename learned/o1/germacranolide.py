"""
Classifies: CHEBI:73011 germacranolide
"""
from rdkit import Chem

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    A germacranolide is a sesquiterpene lactone based on a germacrane skeleton,
    which consists of a macrocyclic ring (typically 10-membered) fused with a lactone ring.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Find macrocycles (rings of size >=9)
    macrocycles = [set(ring) for ring in atom_rings if len(ring) >= 9]
    if not macrocycles:
        return False, "No macrocyclic ring (size >=9) found"
    
    # Identify lactone functional groups (cyclic esters)
    # Lactone pattern: O=C-O in a ring
    lactone_pattern = Chem.MolFromSmarts("C(=O)O")
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    
    # Filter lactone matches that are in a ring
    lactone_in_ring = False
    for match in lactone_matches:
        if all(mol.GetAtomWithIdx(idx).IsInRing() for idx in match):
            lactone_in_ring = True
            break
    if not lactone_in_ring:
        return False, "No lactone functional group found in a ring"
    
    # Optionally, check if lactone is part of macrocycle
    fused = False
    for macrocycle in macrocycles:
        for match in lactone_matches:
            # Check if lactone atoms are within the macrocycle ring
            if set(match).issubset(macrocycle):
                fused = True
                break
        if fused:
            break
    if not fused:
        return False, "Lactone ring not part of macrocyclic structure"
    
    # Optional: Check for alpha-methylene-gamma-lactone moiety
    alpha_methylene_gamma_lactone = Chem.MolFromSmarts("C=CC(=O)O")
    if not mol.HasSubstructMatch(alpha_methylene_gamma_lactone):
        return False, "No alpha-methylene-gamma-lactone moiety found"
    
    return True, "Molecule is a germacranolide (macrocyclic lactone based on germacrane skeleton)"