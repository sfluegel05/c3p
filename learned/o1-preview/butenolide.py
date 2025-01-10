"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: butenolide
"""
from rdkit import Chem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is a gamma-lactone that consists of a 2-furanone skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information
    ri = mol.GetRingInfo()
    atoms = mol.GetAtoms()
    
    # Iterate over all rings
    for ring in ri.BondRings():
        # Check if ring is five-membered
        if len(ring) != 5:
            continue
        
        # Get atoms in the ring
        ring_atoms = set()
        for bond_idx in ring:
            bond = mol.GetBondWithIdx(bond_idx)
            ring_atoms.add(bond.GetBeginAtomIdx())
            ring_atoms.add(bond.GetEndAtomIdx())
        
        # Count oxygens in the ring
        o_in_ring = 0
        for idx in ring_atoms:
            atom = atoms[idx]
            if atom.GetSymbol() == 'O':
                o_in_ring += 1
        if o_in_ring != 1:
            continue  # Proceed if exactly one oxygen in ring
        
        # Check for carbonyl group in ring (C=O)
        carbonyl_in_ring = False
        for idx in ring_atoms:
            atom = atoms[idx]
            if atom.GetSymbol() == 'C':
                for nbr in atom.GetNeighbors():
                    if nbr.GetSymbol() == 'O' and nbr.GetIdx() in ring_atoms:
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            carbonyl_in_ring = True
                            break
        if not carbonyl_in_ring:
            continue  # Proceed if carbonyl group is part of ring
        
        # Check for at least one double bond in ring
        double_bond_in_ring = False
        for bond_idx in ring:
            bond = mol.GetBondWithIdx(bond_idx)
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                double_bond_in_ring = True
                break
        if not double_bond_in_ring:
            continue  # Proceed if ring contains a double bond
        
        # If all conditions are met
        return True, "Contains butenolide ring (gamma-lactone with 2-furanone skeleton)"
    
    # No matching ring found
    return False, "No butenolide ring found"