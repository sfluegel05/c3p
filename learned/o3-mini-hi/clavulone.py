"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: Clavulone – a class of esterified prostanoids obtained from marine corals.
Heuristic criteria (improved from previous attempt):
  1. The molecule must contain a non‐fused five–membered ring (isolated cyclopentenone)
     that has exactly one exocyclic carbonyl group (a C=O bond from a ring atom) 
     and exactly one C=C double bond inside the ring.
  2. At least one ester substituent (–O–C(=O)–R) must be directly attached to a ring atom.
These criteria aim to capture the key scaffold of clavulones while 
reducing false positives from polycyclic or unrelated systems.
"""

from rdkit import Chem

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone (esterified prostanoid) based on its SMILES string.
    
    Heuristic criteria:
      1. The molecule contains a non-fused 5-membered ring (a candidate cyclopentenone)
         with exactly one exocyclic carbonyl group (a double bond from a ring carbon to an oxygen
         that is not in the ring) and exactly one non-carbonyl double bond between ring atoms.
      2. At least one ester group is attached directly to one of the cyclopentenone ring atoms.
         (That is, an oxygen substituent that is further connected to a carbonyl group.)
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a clavulone by our heuristic, False otherwise.
        str : An explanation for the classification decision.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    candidate_ring = None
    
    # Find a 5-membered ring that may correspond to a cyclopentenone scaffold.
    for ring in all_rings:
        if len(ring) != 5:
            continue  # only interested in 5-membered rings
        
        # Check if this ring is "non-fused": every atom appears in only one ring.
        is_non_fused = True
        for atom_idx in ring:
            count = sum([1 for r in all_rings if atom_idx in r])
            if count > 1:
                is_non_fused = False
                break
        if not is_non_fused:
            continue
        
        # Count exocyclic carbonyls:
        # For each atom in the ring, look for a neighbor that is oxygen (atomic num 8)
        # not in the ring and is double-bonded.
        exo_carbony_count = 0
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:  # oxygen
                    bond = mol.GetBondBetweenAtoms(atom_idx, nbr.GetIdx())
                    if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                        exo_carbony_count += 1
        # We require exactly one exocyclic carbonyl from the ring.
        if exo_carbony_count != 1:
            continue
        
        # Count the double bonds internal to the ring (between two ring atoms).
        double_bond_count = 0
        # look at bonds connecting two atoms in the ring (each bond only once)
        ring_atoms_set = set(ring)
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if (a1 in ring_atoms_set and a2 in ring_atoms_set and 
                bond.GetBondType() == Chem.BondType.DOUBLE):
                double_bond_count += 1
        # In a cyclopentenone, besides the carbonyl double bond (which is exocyclic),
        # we expect exactly one additional (carbon-carbon) double bond in the ring.
        if double_bond_count != 1:
            continue
        
        # If we reach here, this ring looks like a cyclopentenone candidate.
        candidate_ring = ring
        break
    
    if candidate_ring is None:
        return False, "No suitable non-fused cyclopentenone ring found (5-membered ring with one exocyclic carbonyl and one internal C=C)"
    
    # Criterion 2: Check for at least one ester substituent attached to the cyclopentenone ring.
    # For each atom in the candidate ring, check its neighbors that lie outside the ring.
    ester_found = False
    for atom_idx in candidate_ring:
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in candidate_ring:
                continue
            # We expect the ester to be attached via an oxygen.
            if nbr.GetAtomicNum() != 8:
                continue
            # Typically this O is singly bonded to the ring.
            bond = mol.GetBondBetweenAtoms(atom_idx, nbr.GetIdx())
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            # Now check if that oxygen is connected (besides the ring atom) to a carbonyl carbon.
            for o_nbr in nbr.GetNeighbors():
                if o_nbr.GetIdx() == atom_idx:
                    continue
                if o_nbr.GetAtomicNum() != 6:
                    continue
                bond2 = mol.GetBondBetweenAtoms(nbr.GetIdx(), o_nbr.GetIdx())
                if bond2 is not None and bond2.GetBondType() == Chem.BondType.DOUBLE:
                    # check that the neighbor of this carbonyl is oxygen (the C=O)
                    for c_nbr in o_nbr.GetNeighbors():
                        if c_nbr.GetIdx() == nbr.GetIdx():
                            continue
                        if c_nbr.GetAtomicNum() == 8:
                            ester_found = True
                            break
                if ester_found:
                    break
            if ester_found:
                break
        if ester_found:
            break
    
    if not ester_found:
        return False, "Cyclopentenone ring found but no ester substituent attached to it."
    
    return True, "Contains a non-fused cyclopentenone ring with an ester substituent – consistent with clavulone"

# Example usage:
if __name__ == "__main__":
    # Example: punaglandin 2 (a valid clavulone)
    smiles_example = "ClC=1C(=O)[C@@]([C@@](O)(C/C=C\\CCCCC)C1)([C@@H](OC(=O)C)[C@H](OC(=O)C)[C@@H](OC(=O)C)CCCC(OC)=O)[H]"
    result, reason = is_clavulone(smiles_example)
    print(result, reason)