"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: CHEBI:50539 butenolide
Butenolide is a gamma-lactone with a 2-furanone skeleton (five-membered ring with one oxygen, a carbonyl group, and conjugated double bonds)
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide contains a five-membered lactone ring (gamma-lactone) with a 2-furanone core structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Get all rings in the molecule
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    for ring in rings:
        if len(ring) != 5:  # Check for five-membered ring
            continue
        
        oxygen_count = 0
        carbonyl_in_ring = False
        
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Count oxygen atoms in the ring
            if atom.GetAtomicNum() == 8:
                oxygen_count += 1
            
            # Check for carbonyl group (C=O) in the ring
            if atom.GetAtomicNum() == 6:
                for bond in atom.GetBonds():
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        neighbor = bond.GetOtherAtom(atom)
                        if neighbor.GetAtomicNum() == 8:
                            carbonyl_in_ring = True
        
        # Check if exactly one oxygen and at least one carbonyl group exists in the ring
        if oxygen_count == 1 and carbonyl_in_ring:
            # Additional check for at least one double bond in the ring (conjugation)
            has_conjugated_double = any(
                bond.GetBondType() == Chem.BondType.DOUBLE
                for bond in mol.GetBonds()
                if bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring
            )
            
            if has_conjugated_double:
                return True, "Contains a five-membered lactone ring with conjugated carbonyl system"
    
    return False, "No five-membered lactone ring with required features found"