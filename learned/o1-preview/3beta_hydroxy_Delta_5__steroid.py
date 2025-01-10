"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
"""
Classifies: 3beta-hydroxy-Delta(5)-steroid

Definition: 'Any 3beta-hydroxy-steroid that contains a double bond between positions 5 and 6.'
"""
from rdkit import Chem

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    A 3beta-hydroxy-Delta(5)-steroid is any 3beta-hydroxy-steroid that contains a double bond between positions 5 and 6 (Delta5).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define steroid nucleus pattern (four fused rings with sizes 6-6-6-5)
    steroid_nucleus = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4C3=CC=C4')  # Simplified steroid nucleus
    if steroid_nucleus is None:
        return False, "Invalid steroid nucleus SMARTS pattern"
    if not mol.HasSubstructMatch(steroid_nucleus):
        return False, "No steroid nucleus found"
    
    # Define 3beta-hydroxy group pattern
    hydroxy_3beta = Chem.MolFromSmarts('[C@H]([C;R0][OH])')  # Chiral center connected to hydroxyl group
    if hydroxy_3beta is None:
        return False, "Invalid hydroxy_3beta SMARTS pattern"
    if not mol.HasSubstructMatch(hydroxy_3beta):
        return False, "No 3beta-hydroxy group found"
    
    # Define Delta(5) double bond between positions 5 and 6
    delta5_double_bond = Chem.MolFromSmarts('C=C')  # Double bond between two carbons
    if delta5_double_bond is None:
        return False, "Invalid delta5_double_bond SMARTS pattern"
    matches = mol.GetSubstructMatches(delta5_double_bond)
    if not matches:
        return False, "No double bonds found"
    # Check if the double bond is between positions 5 and 6 in the steroid nucleus
    # Since atom numbering may not correspond, we approximate by checking if the double bond is in ring C
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    found_double_bond_in_ring = False
    for bond in mol.GetBonds():
        if bond.IsInRing() and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            # Check if both atoms are in the same ring of size 6 (ring C)
            for ring in atom_rings:
                if len(ring) == 6 and begin_atom.GetIdx() in ring and end_atom.GetIdx() in ring:
                    found_double_bond_in_ring = True
                    break
            if found_double_bond_in_ring:
                break
    if not found_double_bond_in_ring:
        return False, "No Delta(5) double bond found between positions 5 and 6"
    
    return True, "Molecule is a 3beta-hydroxy-Delta(5)-steroid"