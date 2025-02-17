"""
Classifies: CHEBI:46640 diketone
"""
"""
Classifies: A compound that contains exactly two ketone functionalities (diketone).
A ketone functionality (R-CO-R) is defined here as a carbon atom that:
  - Is not aromatic (and—if in a ring—it is not part of a fully aromatic ring),
  - Has a double bond to exactly one oxygen,
  - Has no attached hydrogen (so as not to be an aldehyde),
  - And is connected via two single bonds (sigma bonds) to carbon atoms.
This protocol is designed to count only true “ketone” groups and ignore carbonyl groups in 
aromatic (quinone‐like) systems.
"""

from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone (contains exactly two independent ketone functionalities)
    based on its SMILES string.
    
    A ketone functionality is detected by looking for a carbon (non-aromatic) that has exactly one 
    double bond to an oxygen, no attached hydrogen, and two single bonds to carbon atoms. For atoms 
    that are in rings, we further check: if the carbonyl carbon is in any fully aromatic ring then 
    it is ignored.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule has exactly two ketone functionalities; False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ketone_count = 0

    # For ring aromaticity check, get ring information once.
    ring_info = mol.GetRingInfo().AtomRings()
    
    for atom in mol.GetAtoms():
        # We are looking for a carbon atom
        if atom.GetSymbol() != "C":
            continue
        # Exclude if the carbon itself is flagged as aromatic.
        if atom.GetIsAromatic():
            continue
        
        # If the atom is in any ring that is fully aromatic, then skip it.
        atom_idx = atom.GetIdx()
        skip_due_to_aromatic_ring = False
        for ring in ring_info:
            if atom_idx in ring:
                # If all atoms in the ring are aromatic, consider the carbonyl as embedded in an aromatic ring.
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    skip_due_to_aromatic_ring = True
                    break
        if skip_due_to_aromatic_ring:
            continue
        
        # Exclude carbon atoms that have any hydrogen (would be aldehyde if attached to a carbonyl)
        if atom.GetTotalNumHs() != 0:
            continue

        bonds = atom.GetBonds()
        # Identify if this carbon has exactly one double bond to oxygen.
        doubleO_bonds = []
        single_bonds = []
        for bond in bonds:
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() == 8:
                    doubleO_bonds.append(bond)
            elif bond.GetBondType() == Chem.BondType.SINGLE:
                single_bonds.append(bond)
        # Must be exactly one double bond to oxygen.
        if len(doubleO_bonds) != 1:
            continue
        
        # For a (true) ketone the carbonyl carbon should have exactly two single bonds.
        if len(single_bonds) != 2:
            continue
        
        # Both single-bond neighbors must be carbon atoms.
        if any(bond.GetOtherAtom(atom).GetAtomicNum() != 6 for bond in single_bonds):
            continue
        
        # Passed all tests, count this as one ketone functionality.
        ketone_count += 1

    if ketone_count == 2:
        return True, "Contains exactly two ketone functionalities"
    else:
        return False, f"Found {ketone_count} ketone functionalities; exactly 2 are required for a diketone"

# Example usage; you can test using one of the provided diketone SMILES.
if __name__ == "__main__":
    # Example: 2,3-Undecanedione SMILES "O=C(CCCCCCCC)C(=O)C"
    test_smiles = "O=C(CCCCCCCC)C(=O)C"
    result, reason = is_diketone(test_smiles)
    print(f"Result: {result}\nReason: {reason}")