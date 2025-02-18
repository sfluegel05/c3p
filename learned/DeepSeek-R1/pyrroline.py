"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: CHEBI:36437 pyrroline
"""
from rdkit import Chem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is a dihydropyrrole: a 5-membered monocyclic compound with one double bond and at least one nitrogen atom.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a pyrroline, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for exactly one ring and it's 5-membered
    ring_info = mol.GetRingInfo()
    num_rings = len(ring_info.AtomRings())
    if num_rings != 1:
        return False, f"Found {num_rings} rings (must be monocyclic)"
    
    ring_atoms = ring_info.AtomRings()[0]
    ring_size = len(ring_atoms)
    if ring_size != 5:
        return False, f"Ring size is {ring_size} (must be 5-membered)"
    
    # Check for exactly one double bond in the ring
    ring_bonds = []
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in ring_atoms and bond.GetEndAtomIdx() in ring_atoms:
            ring_bonds.append(bond)
    
    double_bonds = sum(1 for bond in ring_bonds if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds != 1:
        return False, f"Found {double_bonds} double bonds in ring (needs 1)"
    
    # Check for at least one nitrogen in the ring
    has_nitrogen = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring_atoms)
    if not has_nitrogen:
        return False, "No nitrogen atom in the ring"
    
    # Ensure it's a heterocycle (no carbon-only rings)
    # Already checked for nitrogen, so this is covered
    
    return True, "5-membered monocyclic ring with one double bond and a nitrogen"