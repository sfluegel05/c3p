"""
Classifies: CHEBI:73011 germacranolide
"""
"""
Classifies: CHEBI:XXXXX germacranolide
"""
from rdkit import Chem

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    A germacranolide is a sesquiterpene lactone derived from the germacrane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Identify gamma-lactone rings (5-membered with O and C=O)
    lactone_rings = []
    for ring in atom_rings:
        if len(ring) != 5:
            continue
        has_oxygen = False
        has_carbonyl = False
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 8:
                has_oxygen = True
            if atom.GetAtomicNum() == 6:
                for bond in atom.GetBonds():
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        neighbor = bond.GetOtherAtom(atom)
                        if neighbor.GetAtomicNum() == 8:
                            has_carbonyl = True
        if has_oxygen and has_carbonyl:
            lactone_rings.append(ring)
    
    if not lactone_rings:
        return False, "No gamma-lactone ring detected"
    
    # Check for fusion between lactone and 10-membered ring
    for lactone_ring in lactone_rings:
        lactone_bonds = set()
        for i in range(len(lactone_ring)):
            a = lactone_ring[i]
            b = lactone_ring[(i+1) % len(lactone_ring)]
            bond = mol.GetBondBetweenAtoms(a, b)
            if bond:
                lactone_bonds.add(bond.GetIdx())
        
        for other_ring in atom_rings:
            if len(other_ring) != 10:
                continue
            other_bonds = set()
            for i in range(len(other_ring)):
                a = other_ring[i]
                b = other_ring[(i+1) % len(other_ring)]
                bond = mol.GetBondBetweenAtoms(a, b)
                if bond:
                    other_bonds.add(bond.GetIdx())
            
            if lactone_bonds & other_bonds:
                return True, "Gamma-lactone fused to 10-membered germacrane skeleton"
    
    return False, "No fused gamma-lactone and 10-membered ring found"