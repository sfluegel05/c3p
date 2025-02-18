"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: isoflavones – any isoflavonoid with a 3-aryl-1-benzopyran-4-one (3-aryl-4H-chromen-4-one) skeleton and its substituted derivatives.
Algorithm:
  1. Parse the molecule from the given SMILES.
  2. Loop over all rings (using RDKit ring info) looking for a six‐membered candidate “pyran” ring that:
       • contains exactly one oxygen (the hetero atom in the pyran ring) 
       • has at least one carbon that is part of a carbonyl (a double‐bonded oxygen, which should be at position 4).
  3. Identify a fused benzene ring. In flavonoids the benzopyran core is fused to an aromatic (benzene) ring. 
     We mark those atoms shared with a fully aromatic six‐membered ring.
  4. In the candidate (“pyran”) ring, atoms that are neither the oxygen nor the carbonyl and are not part of 
     the fusion with the benzene ring are potential “3‑position” atoms.
  5. Check whether any of these candidate atoms carries an external substituent that is itself part of 
     a fully aromatic six‐membered ring (this is taken as the “3‑aryl” group).
  6. If so, return True with a message; otherwise return False.
  
Note: This method is heuristic and may fail in borderline cases.
"""

from rdkit import Chem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone (3-aryl-1-benzopyran-4-one) based on its SMILES string.
    
    Args:
         smiles (str): SMILES string of the molecule.
         
    Returns:
         bool: True if the molecule is classified as an isoflavone, False otherwise.
         str: Reason explaining the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information for the whole molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Loop over rings looking for a candidate six-membered pyran ring.
    for ring in ring_info:
        if len(ring) != 6:
            continue
        # Count the number of oxygen atoms in the ring.
        oxygen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if oxygen_count != 1:
            continue
        
        # Look for a carbon in the ring that is part of a carbonyl (C=O)
        has_carbonyl = False
        carbonyl_atom = None
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                # Check if the neighbor is oxygen and the bond is a double bond.
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond is not None and bond.GetBondType().name == "DOUBLE":
                        has_carbonyl = True
                        carbonyl_atom = atom.GetIdx()  # mark the carbonyl carbon (expected to be the 4-position)
                        break
            if has_carbonyl:
                break
        if not has_carbonyl:
            continue
        
        candidate_ring = set(ring)
        
        # Identify fused benzene rings. In an isoflavone the benzopyran core is fused to a benzene ring.
        fused_atoms = set()
        for other_ring in ring_info:
            if set(other_ring) == candidate_ring:
                continue
            if len(other_ring) != 6:
                continue
            # Check if every atom in the other ring is aromatic.
            if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in other_ring):
                continue
            intersection = candidate_ring.intersection(other_ring)
            if len(intersection) >= 2:
                fused_atoms.update(intersection)
                
        # Identify candidate positions on the pyran ring for the 3-aryl substituent.
        # We do not want the oxygen or the carbonyl carbon, nor atoms in the fusion with the benzene ring.
        candidate_positions = []
        for idx in candidate_ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:  # skip oxygen
                continue
            if carbonyl_atom is not None and idx == carbonyl_atom:
                continue
            if idx in fused_atoms:
                continue
            candidate_positions.append(idx)
        
        # For each candidate position, check if it has an external aromatic substituent.
        for pos in candidate_positions:
            atom = mol.GetAtomWithIdx(pos)
            for nbr in atom.GetNeighbors():
                # Skip neighbor if it belongs to the candidate (pyran) ring.
                if nbr.GetIdx() in candidate_ring:
                    continue
                if not nbr.GetIsAromatic():
                    continue
                # Look if this neighbor is part of a fully aromatic six-membered ring.
                for ring2 in ring_info:
                    if len(ring2) == 6 and nbr.GetIdx() in ring2 and all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring2):
                        return True, "Molecule contains a 3-aryl-1-benzopyran-4-one skeleton and its substituted derivatives."
    return False, "Molecule does not contain a suitable 3-aryl-1-benzopyran-4-one skeleton."

# Example test cases – note that the outcome on borderline examples may vary with heuristic methods.
if __name__ == "__main__":
    test_smiles = [
        # Expected to be isoflavones:
        "COc1ccc(cc1)-c1coc2cc(O)ccc2c1=O",  # formononetin
        "O1C2=C(CC=C(C)C)C(O)=CC(O)=C2C(=O)C(C3=CC=C(O)C=C3)=C1",  # Lupiwighteone
        "COc1cc2c(cc1O)occ(-c1ccc(O)cc1)c2=O",  # glycitein
        # Expected negatives:
        "O=C1OC2=C(C(O)=C(O)C=C2C=3C1=C(O)C=C(OC)C3)C4=C(O)C(O)=CC5=C4C6=C(C(O)=CC(=C6)OC)C(O5)=O",  # Verrulactone B (false positive previously)
    ]
    for smi in test_smiles:
        result, reason = is_isoflavones(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")