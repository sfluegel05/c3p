"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
"""
Classifies: Iridoid monoterpenoid (and secoiridoid derivatives)
Definition:
  Iridoid monoterpenoids are biosynthetically derived from isoprene units and 
  typically possess a fused bicyclic core – a cyclopentane ring fused to a six-membered 
  oxygen heterocycle. In some cases, ring opening (secoiridoids) is observed.
  
  Our algorithm has two stages:
    Stage 1 – Fused bicyclic core detection.
      * Look for one 5-membered ring that is composed only of nonaromatic carbons.
      * Look for one 6-membered ring that contains exactly one (nonaromatic) oxygen 
        and five carbons (all nonaromatic).
      * The two rings must share exactly 2 atoms.
      * Moreover, the two shared atoms must be connected (i.e. they are consecutive in 
        the atom order in both rings) so that the fusion is “real”.
    Stage 2 – Rescue for secoiridoids.
      * Look for any nonaromatic cyclopentane (5-membered, all carbons) that carries at 
        least 2 exocyclic carbonyl substituents.
"""

from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid / secoiridoid derivative
    based on its SMILES string.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): Tuple; whether classified as iridoid monoterpenoid and the reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"

    # Helper function to check if in a ring ordering the two indices are adjacent.
    def are_adjacent_in_ring(ring_tuple, a, b):
        n = len(ring_tuple)
        for i in range(n):
            # Using wrap-around (neighbor indices in ring)
            if (ring_tuple[i] == a and ring_tuple[(i+1)%n] == b) or \
               (ring_tuple[i] == b and ring_tuple[(i+1)%n] == a):
                return True
        return False

    # Stage 1: Check for fused bicyclic core.
    # Iterate over all pairs of rings with sizes 5 and 6.
    for i in range(len(ring_info)):
        for j in range(i+1, len(ring_info)):
            ring1 = ring_info[i]
            ring2 = ring_info[j]
            size1, size2 = len(ring1), len(ring2)
            if sorted([size1, size2]) != [5, 6]:
                continue

            # Identify which ring is 5-membered and which is 6-membered.
            if size1 == 5:
                ring5 = ring1
                order5 = ring1  # preserving order
                ring6 = ring2
                order6 = ring2
            else:
                ring5 = ring2
                order5 = ring2
                ring6 = ring1
                order6 = ring1

            # Check that the 5-membered ring:
            valid_ring5 = True
            for idx in ring5:
                atom = mol.GetAtomWithIdx(idx)
                # Must be carbon, not aromatic. (We allow any sp3 or even non-specified hybridization.)
                if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
                    valid_ring5 = False
                    break
            if not valid_ring5:
                continue

            # Check the 6-membered ring: count oxygen and carbon atoms.
            count_C = 0
            count_O = 0
            valid_ring6 = True
            for idx in ring6:
                atom = mol.GetAtomWithIdx(idx)
                # Disallow aromatic atoms.
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

            # Check for exactly two shared atoms.
            shared_atoms = set(ring5).intersection(set(ring6))
            if len(shared_atoms) != 2:
                continue

            # Check that the two shared atoms are consecutive in both rings.
            shared_atoms = list(shared_atoms)
            if not (are_adjacent_in_ring(order5, shared_atoms[0], shared_atoms[1]) and
                    are_adjacent_in_ring(order6, shared_atoms[0], shared_atoms[1])):
                continue

            reason = ("Found fused bicyclic core: a 5-membered cyclopentane (nonaromatic, all C) and "
                      "a 6-membered oxygen heterocycle (1 O and 5 C, nonaromatic) share exactly 2 atoms, "
                      "with the shared atoms bonded in both rings.")
            return True, reason

    # Stage 2: Rescue for secoiridoids.
    # Look for any 5-membered ring made entirely of nonaromatic carbons.
    for ring in ring_info:
        if len(ring) != 5:
            continue
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and 
                   not mol.GetAtomWithIdx(idx).GetIsAromatic()
                   for idx in ring):
            continue
        
        # For each atom in the ring, check exocyclic bonds that are carbonyls.
        carbonyl_count = 0
        # Use a set to avoid double-counting the same substituent.
        counted_substituents = set()
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                bond = mol.GetBondBetweenAtoms(idx, nbr.GetIdx())
                if bond is None or bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue
                # The neighbor should be a carbon that is doubly bound to oxygen.
                if nbr.GetAtomicNum() != 6:
                    continue
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetIdx() == idx:
                        continue
                    bond2 = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                    if bond2 is not None and bond2.GetBondType() == Chem.BondType.DOUBLE and nbr2.GetAtomicNum() == 8:
                        # Count only if this exocyclic carbon (nbr) has not been already counted.
                        if nbr.GetIdx() not in counted_substituents:
                            carbonyl_count += 1
                            counted_substituents.add(nbr.GetIdx())
                        break  # count each substituent only once
        if carbonyl_count >= 2:
            reason = ("Found cyclopentane ring (all nonaromatic C) with {} exocyclic carbonyl substituents, "
                      "suggestive of a secoiridoid skeleton.".format(carbonyl_count))
            return True, reason

    return False, "No fused iridoid or secoiridoid core detected based on our heuristics."

# Example usage for testing:
if __name__ == "__main__":
    examples = [
        # True positive examples (should be classified as iridoid monoterpenoid)
        ("O1C(OC)C2C([C@H](OC)C[C@@]2(O)C)CC1OC", "Mioporosidegenin"),
        ("O1[C@@H]2C=C(CC[C@]23[C@@H](CC[C@H]3C)[C@](C1)(O)C)C", "Cordycepol A"),
        # Secoiridoid rescue example:
        ("CC1CCC(C(C)C=O)C1C=O", "iridodial"),
        # Examples known to be iridoid monoterpenoids / secoiridoids (false negatives previously):
        ("[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O[C@H]2/C(/[C@](C(=CO2)C(=O)OC)([H])CC(=O)OCCC=3C=CC(=C(C3)O)O)=C/C", "oleuropein"),
        ("[C@]12([C@@]([C@H](CC1)C)(C(OC=C2C)=O)[H])[H]", "cis-trans-nepetalactone"),
    ]
    
    for smi, name in examples:
        flag, msg = is_iridoid_monoterpenoid(smi)
        print(f"{name}: {flag} => {msg}")