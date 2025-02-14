"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (35 <= c_count <= 50):
        return False, f"Carotenoids typically have around 40 carbons, found {c_count}"

    # Look for conjugated double bonds (more flexible count)
    conjugated_pattern = Chem.MolFromSmarts("C=C")
    conjugated_bonds = mol.GetSubstructMatches(conjugated_pattern)
    if len(conjugated_bonds) < 8:
        return False, "Insufficient conjugated double bonds"
    
    # Oxygen atoms presence check
    num_oxygen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygen < 1:
        return False, "No oxygen atoms found"

    # Check the ring presence; allow linear polyene structure
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:  
        return True, "Cyclic structure matches carotenoid properties"
    else:
        # Allow linear systems with conjugated double bonds.
        if len(conjugated_bonds) >= 8 and num_oxygen >= 1:
            return True, "Linear conjugated structure with oxygen matches carotenoid properties"

    return False, "Structure does not match typical carotenoid characteristics"