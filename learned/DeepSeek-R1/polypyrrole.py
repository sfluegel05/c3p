"""
Classifies: CHEBI:38077 polypyrrole
"""
"""
Classifies: CHEBI:xxxxx polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole contains two or more pyrrole units (5-membered aromatic rings with one nitrogen).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypyrrole, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    pyrrole_count = 0
    rings = mol.GetRingInfo().AtomRings()
    
    for ring in rings:
        if len(ring) != 5:
            continue  # Not a 5-membered ring
        # Check if all atoms in the ring are aromatic
        aromatic = True
        for atom_idx in ring:
            if not mol.GetAtomWithIdx(atom_idx).GetIsAromatic():
                aromatic = False
                break
        if not aromatic:
            continue
        # Count nitrogen atoms in the ring
        nitrogen_count = sum(1 for atom_idx in ring if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 7)
        if nitrogen_count == 1:
            pyrrole_count += 1
    
    if pyrrole_count >= 2:
        return True, f"Contains {pyrrole_count} pyrrole rings"
    else:
        return False, f"Found {pyrrole_count} pyrrole rings, need at least 2"