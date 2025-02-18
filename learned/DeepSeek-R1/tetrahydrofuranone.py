"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
"""
Classifies: CHEBI:????? tetrahydrofuranone
"""
from rdkit import Chem
from rdkit.Chem import Mol

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is an oxolane (tetrahydrofuran) with an oxo substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for five-membered rings with exactly one oxygen atom
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) != 5:
            continue
        o_count = sum(1 for atom_idx in ring if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 8)
        if o_count != 1:
            continue
        
        # Check for carbonyl group (C=O) in the ring or as a substituent
        # Case 1: Carbonyl in the ring
        carbonyl_in_ring = False
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:
                continue
            for bond in atom.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    neighbor = bond.GetOtherAtom(atom)
                    if neighbor.GetAtomicNum() == 8:
                        carbonyl_in_ring = True
                        break
            if carbonyl_in_ring:
                break
        if carbonyl_in_ring:
            return True, "Five-membered oxolane ring with carbonyl group"
        
        # Case 2: Carbonyl as a substituent
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    for bond in neighbor.GetBonds():
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            other = bond.GetOtherAtom(neighbor)
                            if other.GetAtomicNum() == 8 and neighbor.GetIdx() not in ring:
                                return True, "Five-membered oxolane ring with carbonyl substituent"
    
    return False, "No five-membered oxolane ring with oxo substituent found"