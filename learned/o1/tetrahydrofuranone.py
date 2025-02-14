"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is defined as 'Any oxolane having an oxo- substituent at any position on the tetrahydrofuran ring.'
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string to create an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    found = False  # Flag to indicate if tetrahydrofuranone ring is found

    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Iterate over all rings in the molecule
    for ring in rings:
        # Consider only five-membered rings
        if len(ring) != 5:
            continue
        
        # Count the number of oxygen atoms in the ring
        oxygen_in_ring = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        if len(oxygen_in_ring) != 1:
            continue  # Not an oxolane ring
        
        # Identify ring carbons that are double-bonded to oxygen (carbonyls)
        carbonyl_carbons = []
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:  # Only carbon atoms
                continue
            for bond in atom.GetBonds():
                if bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() == 8:  # Oxygen
                    carbonyl_carbons.append(atom_idx)
                    break  # Found a carbonyl, move to next atom
        if len(carbonyl_carbons) == 1:
            # Exactly one carbonyl in the ring, classify as tetrahydrofuranone
            found = True
            break

    if found:
        return True, "Contains an oxolane ring with a carbonyl group (tetrahydrofuranone)"
    else:
        return False, "Does not contain an oxolane ring with a carbonyl group"