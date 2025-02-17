"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
"""
Classifies: Iridoid monoterpenoid
Definition:
  Iridoid monoterpenoids are biosynthesized from isoprene units and typically
  possess a fused bicyclic core – a cyclopentane ring fused to a six-membered
  oxygen heterocycle (usually a tetrahydropyran). In some cases, opening of the
  cyclopentane ring (secoiridoids) is observed. Our improved algorithm uses two stages:
    Stage 1 – Fused bicyclic core detection. We look for a pair of rings:
              • one 5-membered ring composed solely of (nonaromatic, sp3) carbons
              • one 6-membered ring that contains exactly one oxygen (nonaromatic) 
                and five carbons (again sp3)
              • the two rings must share exactly 2 atoms.
    Stage 2 – Rescue for secoiridoids. We look for any cyclopentane ring
              (again, nonaromatic, sp3 carbon only) that carries at least 2 
              exocyclic carbonyl (C=O) substituents.
"""

from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid (or secoiridoid derivative)
    based on its SMILES string.

    The algorithm has two stages:
      Stage 1:
        Look for a fused bicyclic core. For each pair of rings (as returned by RDKit's
        ring detection) where one ring has 5 atoms and the other 6, apply these filters:
         - The 5-membered ring must be composed solely of nonaromatic, sp3 carbons.
         - The 6-membered ring must contain exactly one nonaromatic oxygen and five nonaromatic carbons,
           and those atoms should be sp3 as well.
         - The rings must share exactly 2 atoms.
      Stage 2:
        If no fused core is found, scan every 5-membered nonaromatic cyclopentane (all sp3 carbons) and for 
        each ring count exocyclic carbonyl groups (a double bond from a ring carbon to an exocyclic carbon,
        which in turn is double-bonded to oxygen). If at least 2 such exocyclic carbonyls are attached, 
        classify as a secoiridoid.

    Args:
      smiles (str): SMILES string of the molecule.

    Returns:
      (bool, str): A tuple where the first element is True if the molecule is classified as an
                   iridoid monoterpenoid (or secoiridoid derivative), False otherwise.
                   The second element gives the reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"

    # Stage 1: Look for fused bicyclic core satisfying our improved conditions.
    for i in range(len(ring_info)):
        for j in range(i+1, len(ring_info)):
            ring1 = set(ring_info[i])
            ring2 = set(ring_info[j])
            if (len(ring1) == 5 and len(ring2) == 6) or (len(ring1) == 6 and len(ring2) == 5):
                ring5 = ring1 if len(ring1) == 5 else ring2
                ring6 = ring2 if len(ring2) == 6 else ring1

                # Check that every atom in the 5-membered ring is a nonaromatic carbon with sp3 hybridization.
                valid_ring5 = True
                for idx in ring5:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetAtomicNum() != 6 or atom.GetIsAromatic() or atom.GetHybridization() != Chem.HybridizationType.SP3:
                        valid_ring5 = False
                        break
                if not valid_ring5:
                    continue

                # Check the 6-membered ring: exactly one oxygen and five carbons.
                count_C = 0
                count_O = 0
                valid_ring6 = True
                for idx in ring6:
                    atom = mol.GetAtomWithIdx(idx)
                    # Do not allow aromatic atoms in our ideal core.
                    if atom.GetIsAromatic():
                        valid_ring6 = False
                        break
                    if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.HybridizationType.SP3:
                        count_C += 1
                    elif atom.GetAtomicNum() == 8:
                        # allow oxygen if nonaromatic (hybridization flexibility allowed)
                        count_O += 1
                    else:
                        valid_ring6 = False
                        break
                if not valid_ring6 or count_C != 5 or count_O != 1:
                    continue

                # Check that the rings are fused with exactly 2 shared atoms.
                shared_atoms = ring5.intersection(ring6)
                if len(shared_atoms) == 2:
                    reason = ("Found fused bicyclic core: a 5-membered cyclopentane (nonaromatic, sp3 carbons) and "
                              "a 6-membered ring with one oxygen (5C+1O, nonaromatic) share exactly 2 atoms.")
                    return True, reason

    # Stage 2: Rescue logic for secoiridoids.
    # Look for any 5-membered ring (cyclopentane) that is all carbon, nonaromatic and sp3.
    for ring in ring_info:
        if len(ring) != 5:
            continue
        if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and 
                   not mol.GetAtomWithIdx(idx).GetIsAromatic() and 
                   mol.GetAtomWithIdx(idx).GetHybridization() == Chem.HybridizationType.SP3 
                   for idx in ring):
            continue
        carbonyl_count = 0
        # For each atom in the ring, examine exocyclic bonds.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                bond = mol.GetBondBetweenAtoms(idx, nbr.GetIdx())
                if bond is None or bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue
                # The neighbor should be a carbon (exocyclic carbon) with at least one double bond to oxygen.
                if nbr.GetAtomicNum() != 6:
                    continue
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetIdx() == idx:
                        continue
                    bond2 = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                    if bond2 is not None and bond2.GetBondType() == Chem.BondType.DOUBLE and nbr2.GetAtomicNum() == 8:
                        carbonyl_count += 1
                        break  # count each substituent only once
        if carbonyl_count >= 2:
            reason = ("Found cyclopentane ring (all nonaromatic sp3 C) with {} exocyclic carbonyl substituents "
                      "suggestive of a secoiridoid skeleton.".format(carbonyl_count))
            return True, reason

    return False, "No fused iridoid or secoiridoid core detected based on our heuristics."

# Example usage for testing purposes:
if __name__ == "__main__":
    examples = [
        # True positives based on fused bicyclic core:
        ("[C@]12([C@@]([C@H](CC1)C)(C(OC=C2C)=O)[H])[H]", "cis-trans-nepetalactone"),
        ("O[C@@]12[C@@]([C@@](O)([C@@H](OC(=O)/C=C/C3=CC=C(O)c(O)C=C3)C1)C)(C(OC=C2C(OC)=O)O[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)CO)[H]",
         "Methyl (1s,4ar,6s,7r,7as)-..."),
        ("[H][C@@]12C=CO[C@@H](O[C@]3([H])O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@]1([H])C(COC(=O)C1=CC=C(O)C=C1)=C[C@H]2O",
         "agnuside"),
        # True secoiridoid rescue example:
        ("CC1CCC(C(C)C=O)C1C=O", "iridodial"),
        # A known challenging one (typically false negative with earlier heuristics)
        ("O[C@H]1/C(/[C@](C(=CO1)C(=O)OC)([H])CC(=O)OCCC=2C=CC(=C(C2)O)O)=C/C", "oleuropein aglycone"),
        # A false positive (should be rejected)
        ("O=C1OC2=C([C@](O)(C)CC2=O)C=3C1=C(O)C=C(OC)C3", "Rubralactone"),
    ]
    
    for smi, name in examples:
        flag, msg = is_iridoid_monoterpenoid(smi)
        print(f"{name}: {flag} => {msg}")