"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: Clavulone – a class of esterified prostanoids obtained from marine corals.
Heuristic criteria:
  1. The molecule must contain a 5–membered (cyclopentenone-type) ring that has:
       - exactly one exocyclic carbonyl group (an oxygen outside the ring connected via a double bond to a ring carbon)
       - exactly one double bond between ring atoms.
  2. At least one ester substituent (–O–C(=O)–R) must be directly attached to one of the ring atoms.
These criteria attempt to capture the clavulone scaffold while allowing (fused) rings.
"""

from rdkit import Chem

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone (esterified prostanoid) based on its SMILES string.
    
    Heuristic criteria:
      1. The molecule contains a 5-membered ring (a candidate cyclopentenone) with exactly:
           - one exocyclic carbonyl: a ring atom that is double-bonded to an oxygen that is not part of the ring,
           - one internal (ring–ring) double bond.
      2. At least one ester substituent (an oxygen singly bonded to a ring atom that in turn binds 
         to a carbonyl carbon via a single bond) is found attached to any atom of that ring.
    
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

    # Loop over all 5-membered rings regardless of fusion
    for ring in all_rings:
        if len(ring) != 5:
            continue

        # Count exocyclic carbonyls.
        exo_carbony_count = 0
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Looking for oxygen (atomic num 8) attached by a double bond.
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(atom_idx, nbr.GetIdx())
                    if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                        exo_carbony_count += 1
        
        # We require exactly one exocyclic carbonyl from the ring.
        if exo_carbony_count != 1:
            continue
        
        # Count the internal (ring–ring) double bonds.
        double_bond_count = 0
        ring_atom_set = set(ring)
        # Loop bonds in molecule that connect two atoms in this ring.
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in ring_atom_set and a2 in ring_atom_set:
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    double_bond_count += 1
        
        if double_bond_count != 1:
            continue
        
        # Passed ring criteria so take this as our candidate ring.
        candidate_ring = ring
        break

    if candidate_ring is None:
        return False, ("No suitable cyclopentenone-type ring found: "
                      "need a 5-membered ring with exactly one exocyclic carbonyl "
                      "and one internal double bond.")
    
    # Criterion 2: Look for at least one ester substituent attached to the candidate ring.
    ester_found = False
    for atom_idx in candidate_ring:
        ring_atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in ring_atom.GetNeighbors():
            if nbr.GetIdx() in candidate_ring:
                continue
            # Check that neighbor is oxygen (ester oxygen expected to be outside the ring)
            if nbr.GetAtomicNum() != 8:
                continue
            bond = mol.GetBondBetweenAtoms(atom_idx, nbr.GetIdx())
            # Expect a single bond between ring atom and the oxygen
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            # Now check whether this oxygen is attached (besides the ring) to a carbon that bears a carbonyl.
            for sub_nbr in nbr.GetNeighbors():
                if sub_nbr.GetIdx() == atom_idx:
                    continue
                if sub_nbr.GetAtomicNum() != 6:  # expecting a carbon atom
                    continue
                bond2 = mol.GetBondBetweenAtoms(nbr.GetIdx(), sub_nbr.GetIdx())
                if bond2 is None or bond2.GetBondType() != Chem.BondType.SINGLE:
                    continue
                # Check if this carbon is double-bonded to an oxygen.
                has_carbonyl = False
                for c_nbr in sub_nbr.GetNeighbors():
                    if c_nbr.GetIdx() == nbr.GetIdx():
                        continue
                    if c_nbr.GetAtomicNum() == 8:
                        bond3 = mol.GetBondBetweenAtoms(sub_nbr.GetIdx(), c_nbr.GetIdx())
                        if bond3 is not None and bond3.GetBondType() == Chem.BondType.DOUBLE:
                            has_carbonyl = True
                            break
                if has_carbonyl:
                    ester_found = True
                    break
            if ester_found:
                break
        if ester_found:
            break

    if not ester_found:
        return False, "Cyclopentenone-type ring found but no ester substituent attached to it."
    
    return True, "Contains a cyclopentenone-type 5-membered ring with one exocyclic carbonyl and at least one ester substituent – consistent with clavulone"

# Example usage:
if __name__ == "__main__":
    # Test one of the provided SMILES: punaglandin 2
    smiles_example = "ClC=1C(=O)[C@@]([C@@](O)(C/C=C\\CCCCC)C1)([C@@H](OC(=O)C)[C@H](OC(=O)C)[C@@H](OC(=O)C)CCCC(OC)=O)[H]"
    result, reason = is_clavulone(smiles_example)
    print(result, reason)