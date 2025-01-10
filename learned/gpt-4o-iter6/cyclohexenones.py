"""
Classifies: CHEBI:48953 cyclohexenones
"""
from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is defined as any six-membered alicyclic ketone having one double bond in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Revised SMARTS pattern to match cyclohexenone structure
    # - Ensure six-membered ring
    # - Contains exactly one double bond and one ketone group
    # - Allows for substituents placement
    # - Specific position for the carbonyl group
    cyclohexenone_pattern = Chem.MolFromSmarts("[C,R2]=[C,R2][C,R2][C,R2]C(=O)[C,R2]")
    
    # Check for the cyclohexenone pattern
    if not mol.HasSubstructMatch(cyclohexenone_pattern):
        return False, "No cyclohexenone structure found"
    
    # Verify a single six-membered ring with a double bond
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 6:
            atoms_in_ring = [mol.GetAtomWithIdx(i) for i in ring]
            # Ensure there's one double bond in the ring
            db_count = sum([1 for a in atoms_in_ring if any(nb.GetBondType() == Chem.BondType.DOUBLE for nb in a.GetNeighbors())])
            # Ensure exactly one ketone in the ring
            ketone_count = sum([1 for a in atoms_in_ring if a.GetAtomicNum() == 6 and any(neigh.GetAtomicNum() == 8 for neigh in a.GetNeighbors())])
            if db_count == 1 and ketone_count == 1:
                return True, "Cyclohexenone structure identified with appropriate ring and ketone"
    
    return False, "No six-membered alicyclic ketone with one double bond found in the ring"