"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
"""
Classifies: Iridoid monoterpenoid / secoiridoid derivatives

Definition:
  Iridoid monoterpenoids are derived from isoprene and typically contain a fused bicyclic core.
  The classic core is a cyclopentane ring (nonaromatic, all carbons) fused to a 6‐membered oxygen
  heterocycle (five carbons and one oxygen, nonaromatic) via two shared (and consecutive) carbon atoms.
  Alternatively, secoiridoids are thought to arise by cleavage of a bond in the cyclopentane and are
  marked by a cyclopentane that bears at least two exocyclic carbonyl substituents.

Our revised algorithm has two stages:
  Stage 1 – Look for a fused bicyclic core with a cyclopentane (5-membered all-C ring) and a 6-membered
            ring (5 C and 1 O) that share exactly two atoms (which must be adjacent in both rings). We also
            require that the union of the two rings is made of exactly 8 carbon atoms and 1 oxygen atom.
            
  Stage 2 – If stage 1 fails, “rescue” for secoiridoids by checking for a nonaromatic cyclopentane ring that
            carries at least 2 exocyclic carbonyl substituents. In our updated version, instead of requiring a
            double bond from the ring to the substituent, we check whether the substituent (a carbon) itself
            carries a C=O (a double bond from that carbon to an oxygen).

Note that these heuristics remain simplified and some unusual cases may escape detection.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid / secoiridoid derivative
    based on its SMILES string.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): Tuple of classification result and reason message.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"

    # Helper: check if two atoms (by indices) are adjacent in a given ring ordering.
    def are_adjacent_in_ring(ring_tuple, a, b):
        n = len(ring_tuple)
        for i in range(n):
            if (ring_tuple[i] == a and ring_tuple[(i+1)%n] == b) or \
               (ring_tuple[i] == b and ring_tuple[(i+1)%n] == a):
                return True
        return False

    # Stage 1: Look for fused bicyclic core.
    # Iterate over all pairs of rings. We want one of size 5 and one of size 6.
    for i in range(len(ring_info)):
        for j in range(i+1, len(ring_info)):
            ring1 = ring_info[i]
            ring2 = ring_info[j]
            size1, size2 = len(ring1), len(ring2)
            if sorted([size1, size2]) != [5, 6]:
                continue

            # Identify which ring is 5-membered and which is 6-membered.
            ring5 = ring1 if size1 == 5 else ring2
            order5 = ring5  # ordering as given from ring info
            ring6 = ring1 if size1 == 6 else ring2
            order6 = ring6

            # Check the 5-membered ring: each atom must be a carbon and nonaromatic.
            valid_ring5 = True
            for idx in ring5:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
                    valid_ring5 = False
                    break
            if not valid_ring5:
                continue

            # Check the 6-membered ring: must have exactly 6 atoms and (exactly 5 carbons and 1 oxygen).
            count_C = 0
            count_O = 0
            valid_ring6 = True
            for idx in ring6:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetIsAromatic():
                    valid_ring6 = False
                    break
                if atom.GetAtomicNum() == 6:
                    count_C += 1
                elif atom.GetAtomicNum() == 8:
                    count_O += 1
                else:
                    valid_ring6 = False
                    break
            if not valid_ring6 or count_C != 5 or count_O != 1:
                continue

            # Check that the two rings share exactly two atoms.
            shared_atoms = set(ring5).intersection(set(ring6))
            if len(shared_atoms) != 2:
                continue

            # Verify that the two shared atoms are consecutive (adjacent) in both rings.
            shared_atoms_list = list(shared_atoms)
            if not (are_adjacent_in_ring(order5, shared_atoms_list[0], shared_atoms_list[1]) and
                    are_adjacent_in_ring(order6, shared_atoms_list[0], shared_atoms_list[1])):
                continue

            # Extra check: the union of the fused rings should have exactly 5 (ring5) + 6 (ring6) - 2(shared) = 9 atoms.
            # In the fused core the 5-membered ring gives 5 carbons and the 6-membered ring has 5 carbons and 1 oxygen.
            fused_atoms = set(ring5).union(set(ring6))
            core_C = 0
            core_O = 0
            for idx in fused_atoms:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:
                    core_C += 1
                elif atom.GetAtomicNum() == 8:
                    core_O += 1
            # We expect 8 carbons and 1 oxygen.
            if core_C != 8 or core_O != 1:
                continue
            
            reason = ("Found fused bicyclic core: a 5-membered cyclopentane (nonaromatic, all C) and "
                      "a 6-membered oxygen heterocycle (1 O and 5 C, nonaromatic) share exactly 2 atoms "
                      "that are adjacent in both rings; fused core has 8 C and 1 O.")
            return True, reason
        
    # Stage 2: Rescue for secoiridoids.
    # Look for any 5-membered ring composed solely of nonaromatic carbons.
    for ring in ring_info:
        if len(ring) != 5:
            continue
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and not mol.GetAtomWithIdx(idx).GetIsAromatic() 
                   for idx in ring):
            continue
        
        # For each atom in the ring, check exocyclic substituents.
        # Instead of checking for a double bond in the bond connecting the substituent,
        # we check whether a neighboring atom (not in the ring) is a carbon that itself has
        # a carbonyl (i.e. a double bond to oxygen).
        carbonyl_count = 0
        counted_substituents = set()
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # skip neighbors that are part of the ring
                # To avoid double-counting, use the neighbor index.
                if nbr.GetIdx() in counted_substituents:
                    continue
                # Check that the substituent is a carbon.
                if nbr.GetAtomicNum() != 6:
                    continue
                # Now check if this neighbor carbon is part of a carbonyl group.
                found_carbonyl = False
                for nbr2 in nbr.GetNeighbors():
                    # Skip if nbr2 is the ring atom.
                    if nbr2.GetIdx() == idx:
                        continue
                    bond2 = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                    if bond2 and bond2.GetBondType() == Chem.BondType.DOUBLE and nbr2.GetAtomicNum() == 8:
                        found_carbonyl = True
                        break
                if found_carbonyl:
                    carbonyl_count += 1
                    counted_substituents.add(nbr.GetIdx())
        if carbonyl_count >= 2:
            reason = ("Found cyclopentane ring (nonaromatic, all C) with {} exocyclic carbonyl substituents, "
                      "suggestive of a secoiridoid skeleton.".format(carbonyl_count))
            return True, reason

    return False, "No fused iridoid or secoiridoid core detected based on our heuristics."

# Example usage for testing:
if __name__ == "__main__":
    examples = [
        # True positives (should be classified as iridoid monoterpenoid or secoiridoid)
        ("[C@]12([C@@]([C@H](CC1)C)(C(OC=C2C)=O)[H])[H]", "cis-trans-nepetalactone"),
        ("O[C@@]12[C@@]([C@@](O)([C@@H](OC(=O)/C=C/C3=CC=C(O)c(O)C=C3)C1)C)(C(OC=C2C(OC)=O)O[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)CO)[H]", 
         "Methyl (1s,4ar,6s,7r,7as)-1-(beta-d-glucopyranosyloxy)-4a,7-dihydroxy-6-{[(2e)-3-(4-hydroxyphenyl)-2-propenoyl]oxy}-7-methyl-1,4a,5,6,7,7a-hexahydrocyclopenta[c]pyran-4-carboxylate"),
        ("[H][C@@]12C[C@H](O)[C@H](C)[C@@]1([H])[C@@H](OC=C2C(=O)OC)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O", "loganin"),
        # Secoiridoid rescue example:
        ("CC1CCC(C(C)C=O)C1C=O", "iridodial"),
        # Examples that previously were missed:
        ("[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O[C@H]2/C(/[C@](C(=CO2)C(=O)OC)([H])CC(=O)OCCC=3C=CC(=C(C3)O)O)=C/C", "oleuropein"),
        ("O[C@H]1/C(/[C@](C(=CO1)C(=O)OC)([H])CC(=O)OCCC=2C=CC(=C(C2)O)O)=C/C", "oleuropein aglycone"),
    ]
    
    for smi, name in examples:
        flag, msg = is_iridoid_monoterpenoid(smi)
        print(f"{name}: {flag} => {msg}")