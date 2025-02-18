"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
"""
Classifies: Iridoid monoterpenoid
Definition: Iridoid monoterpenoids are monoterpenoids biosynthesized from isoprene units.
They typically contain a cyclopentane ring fused to a six-membered oxygen heterocycle.
Cleavage of a bond in the cyclopentane ring gives rise to secoiridoids.
This program uses two strategies:
  Stage 1 (ideal iridoid): Look for a fused bicyclic core where one ring is a cyclopentane
         (5 carbons) and the other is a six-membered ring containing exactly one oxygen and 5 carbons.
         In addition, we require that the rings share exactly 2 atoms.
  Stage 2 (rescue logic for secoiridoids): If no fused core is found, then look for a cyclopentane ring
         comprised solely of carbons that carries at least two exocyclic carbonyl (C=O) groups.
Note: These heuristics aim to improve both sensitivity and selectivity.
"""

from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid (or secoiridoid derivative)
    based on its SMILES string.

    Approach:
      Stage 1: Look for a fused bicyclic core. We iterate over each pair of rings (from RDKit's
               ring detection) looking specifically for one 5-membered ring whose atoms are all carbons,
               and one 6-membered ring containing exactly one oxygen (and five carbons).
               In addition, we require that the two rings share exactly 2 atoms.
      Stage 2: If no ideal fused core is found, try to rescue a secoiridoid classification.
               Search for any cyclopentane ring composed solely of carbons. For each atom in the ring,
               count exocyclic substituents that are part of a carbonyl group.
               We classify as secoiridoid if at least 2 such exocyclic carbonyl groups are attached.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an iridoid monoterpenoid (or secoiridoid derivative),
              False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"
    
    # Stage 1: Check for fused iridoid core.
    # Look for two rings, one of size 5 and one of size 6.
    # For the 5-membered ring, all atoms must be carbon.
    # For the 6-membered ring, exactly one atom must be oxygen and the rest carbon.
    # Also, the rings must share exactly 2 atoms.
    for i in range(len(ring_info)):
        for j in range(i+1, len(ring_info)):
            ring1 = set(ring_info[i])
            ring2 = set(ring_info[j])
            # Check for sizes: one must be 5 and the other 6.
            if (len(ring1) == 5 and len(ring2) == 6) or (len(ring1) == 6 and len(ring2) == 5):
                # Identify which is 5-membered and which is 6-membered.
                if len(ring1) == 5:
                    ring5 = ring1
                    ring6 = ring2
                else:
                    ring5 = ring2
                    ring6 = ring1
                
                # Verify the 5-membered ring is all carbon.
                if any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in ring5):
                    continue
                # Verify the 6-membered ring: count carbons and oxygens.
                nC = 0
                nO = 0
                for idx in ring6:
                    atomic_num = mol.GetAtomWithIdx(idx).GetAtomicNum()
                    if atomic_num == 6:
                        nC += 1
                    elif atomic_num == 8:
                        nO += 1
                    else:
                        nC = -1  # if other heteroatoms are present, skip
                        break
                if nC != 5 or nO != 1:
                    continue
                
                # Check shared atoms: require exactly 2 shared atoms.
                shared_atoms = ring5.intersection(ring6)
                if len(shared_atoms) == 2:
                    reason = ("Found fused bicyclic core: a 5-membered cyclopentane (all C) and a 6-membered ring "
                              "with one oxygen (5C+1O) share exactly 2 atoms.")
                    return True, reason
    
    # Stage 2: Rescue approach for secoiridoids.
    # Look for a cyclopentane ring (5-membered ring composed solely of carbons).
    for ring in ring_info:
        if len(ring) == 5:
            if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                continue
            carbonyl_count = 0
            # For each ring atom, examine neighbors not in the ring.
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                for nbr in atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx in ring:
                        continue
                    # For the bond between atom and neighbor, if that neighbor (an exocyclic carbon)
                    # has a double bond to oxygen, count it as a carbonyl.
                    bond = mol.GetBondBetweenAtoms(idx, nbr_idx)
                    if bond is None or bond.GetBondType() != Chem.BondType.DOUBLE:
                        continue
                    # Ensure the neighbor is carbon.
                    if nbr.GetAtomicNum() != 6:
                        continue
                    # Check if this neighbor carbon is in turn double-bonded to at least one oxygen.
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetIdx() == idx:
                            continue
                        bond2 = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                        if bond2 is not None and bond2.GetBondType() == Chem.BondType.DOUBLE and nbr2.GetAtomicNum() == 8:
                            carbonyl_count += 1
                            break  # Count each exocyclic carbonyl only once
            if carbonyl_count >= 2:
                reason = ("Found cyclopentane ring (all C) with {} exocyclic carbonyl substituents "
                          "suggestive of a secoiridoid skeleton.".format(carbonyl_count))
                return True, reason

    return False, "No fused iridoid or secoiridoid core detected based on our heuristics."

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Some example molecules (names given for clarity)
    examples = [
        # True positive based on fused core: cis-trans-nepetalactone
        ("[C@]12([C@@]([C@H](CC1)C)(C(OC=C2C)=O)[H])[H]", "cis-trans-nepetalactone"),
        # A secoiridoid example: iridodial (previously missed)
        ("CC1CCC(C(C)C=O)C1C=O", "iridodial"),
        # A false positive example from prior run:
        ("C1[C@@H](O[C@@H]([C@H]2[C@@H]1C3=C(O2)C=CC(=C3)NS(=O)(=O)C4=CC=CC=C4)CO)CC(=O)NCC(F)(F)F", "Kadsurenone mimic"),
        # A false negative from before: oleuropein aglycone (more challenging).
        ("O[C@H]1/C(/[C@](C(=CO1)C(=O)OC)([H])CC(=O)OCCC=2C=CC(=C(C2)O)O)=C/C", "oleuropein aglycone"),
    ]
    
    for smi, name in examples:
        flag, msg = is_iridoid_monoterpenoid(smi)
        print(f"{name}: {flag} => {msg}")