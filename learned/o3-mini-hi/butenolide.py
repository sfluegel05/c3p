"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: Butenolide – a gamma-lactone consisting of a 2-furanone skeleton and substituted derivatives.
"""
from rdkit import Chem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is defined as a gamma-lactone that consists of a 2-furanone skeleton and its substituted derivatives.
    This means the molecule must contain a five-membered ring (the lactone ring) in which:
      - One ring carbon is double-bonded to an exocyclic oxygen (the lactone carbonyl).
      - There is at least one double bond (unsaturation) among the ring bonds to recall the 2-furanone skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a butenolide, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Obtain ring information from the molecule
    ring_info = mol.GetRingInfo()
    
    # Loop over all rings in the molecule
    for ring in ring_info.AtomRings():
        # We only consider five-membered rings (the gamma-lactone motif must be five-membered)
        if len(ring) != 5:
            continue
        
        lactone_found = False  # Will be True if a ring carbon has an exocyclic C=O (carbonyl)
        unsaturation_found = False  # Will be True if there is at least one double bond within the ring bonds
        
        # Check each atom in the ring for a lactone carbonyl: a carbon atom (atomic number 6)
        # that is double-bonded to an oxygen atom that is not in the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                for nb in atom.GetNeighbors():
                    nb_idx = nb.GetIdx()
                    # Ensure the neighbor is outside the ring and is an oxygen (atomic number 8)
                    if nb_idx not in ring and nb.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb_idx)
                        if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                            lactone_found = True
                            break
            if lactone_found:
                break
        
        # Check for unsaturation (double bond) between adjacent atoms in the ring.
        ring_atoms = list(ring)
        n = len(ring_atoms)
        for i in range(n):
            a1 = ring_atoms[i]
            a2 = ring_atoms[(i+1) % n]  # loop closure between last and first atom
            bond = mol.GetBondBetweenAtoms(a1, a2)
            if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                unsaturation_found = True
                break
        
        # If both criteria are met in a given 5-membered ring, we classify it as a butenolide.
        if lactone_found and unsaturation_found:
            return True, "Contains a 5-membered gamma-lactone ring with a carbonyl and unsaturation (2-furanone core)"
    
    return False, "No suitable 5-membered gamma-lactone (2-furanone) ring detected"