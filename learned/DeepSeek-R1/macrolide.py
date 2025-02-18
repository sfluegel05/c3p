"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolFromSmarts
from rdkit.Chem.rdchem import Mol

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide is a macrocyclic lactone with a ring of twelve or more members.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find oxygen atoms in rings that are part of ester groups (O-C=O)
    ester_pattern = Chem.MolFromSmarts('[O;R]C(=O)')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups in rings"
    
    max_ring_size = 0
    ring_info = mol.GetRingInfo()
    
    for match in ester_matches:
        o_idx, c_idx = match[0], match[1]
        # Check if oxygen and carbonyl carbon are in the same ring
        if not ring_info.AreAtomsInSameRing(o_idx, c_idx):
            continue
        
        # Find all rings containing both atoms and track the largest
        for ring in ring_info.AtomRings():
            if o_idx in ring and c_idx in ring:
                current_size = len(ring)
                if current_size > max_ring_size:
                    max_ring_size = current_size
                break  # Only need the first matching ring for this pair
    
    if max_ring_size >= 12:
        return True, f"Contains {max_ring_size}-membered macrocyclic lactone"
    else:
        return False, f"Largest lactone ring: {max_ring_size} members (needs ≥12)"