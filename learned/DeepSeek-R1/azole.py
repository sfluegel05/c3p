"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: CHEBI:46640 azole
"""
from rdkit import Chem

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.
    An azole is a monocyclic five-membered ring containing at least one nitrogen atom.
    The ring can also include other heteroatoms (N, S, O) and must not be fused with other rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an azole, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Sanitize to ensure proper ring and atom properties
    Chem.SanitizeMol(mol)
    
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    for ring in rings:
        if len(ring) != 5:
            continue  # Not a five-membered ring
        
        # Check for at least one nitrogen in the ring
        has_nitrogen = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring)
        if not has_nitrogen:
            continue
        
        # Check if the ring is monocyclic (none of its atoms are in other rings)
        is_monocyclic = all(ring_info.NumAtomRings(idx) == 1 for idx in ring)
        if is_monocyclic:
            return True, "Contains a five-membered monocyclic ring with at least one nitrogen"
    
    return False, "No five-membered monocyclic ring with nitrogen found"