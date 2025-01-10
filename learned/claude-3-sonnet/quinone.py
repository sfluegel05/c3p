"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: CHEBI:26421 quinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for cyclic ketone pattern (C(=O) in ring)
    ketone_pattern = Chem.MolFromSmarts("[C,c;R]([#6,#1;A])=O")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    if len(ketone_matches) < 2:
        return False, "Must have at least 2 cyclic ketone groups"
    
    # Check if ketones are in the same ring system
    ring_info = mol.GetRingInfo()
    
    # Get ring atoms for each ketone
    ketone_ring_atoms = []
    for match in ketone_matches:
        ketone_carbon = match[0]
        rings_for_ketone = ring_info.AtomRings()
        for ring in rings_for_ketone:
            if ketone_carbon in ring:
                ketone_ring_atoms.append(set(ring))
                break
                
    if not ketone_ring_atoms:
        return False, "Ketone groups must be part of ring system"
    
    # Check if at least 2 ketones share a ring system
    connected_ketones = False
    for i in range(len(ketone_ring_atoms)):
        for j in range(i+1, len(ketone_ring_atoms)):
            if ketone_ring_atoms[i].intersection(ketone_ring_atoms[j]):
                connected_ketones = True
                break
    
    if not connected_ketones:
        return False, "Ketone groups must be in connected ring system"
    
    # Check for conjugation between ketones
    conjugated_pattern = Chem.MolFromSmarts("[C,c;R]([#6,#1;A])=O.[C,c;R]([#6,#1;A])=O")
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "Ketone groups must be conjugated"
    
    # Look for aromatic or conjugated system
    aromatic_atoms = [atom.GetIsAromatic() for atom in mol.GetAtoms()]
    conjugated_bonds = [bond.GetIsConjugated() for bond in mol.GetBonds()]
    
    if not (any(aromatic_atoms) or any(conjugated_bonds)):
        return False, "Must have aromatic or conjugated system"
        
    # Additional check for quinone pattern - alternating single/double bonds
    quinone_pattern = Chem.MolFromSmarts("[C,c;R]([#6,#1;A])=O-[C,c;R]=,:[C,c;R]-[C,c;R]([#6,#1;A])=O")
    if not mol.HasSubstructMatch(quinone_pattern):
        return False, "Must have characteristic quinone pattern"

    return True, "Contains cyclic dione structure with conjugation characteristic of quinones"