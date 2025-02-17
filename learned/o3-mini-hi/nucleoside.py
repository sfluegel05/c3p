"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: Nucleoside
Definition: An N-glycosyl compound that has both a nucleobase 
           (normally adenine, guanine, xanthine, thymine, cytosine or uracil)
           and either a ribose or deoxyribose as functional parents.
Improvement notes:
  - Locate a candidate sugar ring as a non-aromatic 5-membered ring with 4 carbons and 1 oxygen.
  - Identify a bond from a sugar atom (ideally the anomeric carbon) to a nitrogen atom that is 
    not part of the sugar but belongs to a ring of size 5 or 6 as a minimal nucleobase.
  - Relax full aromaticity for the base so that partially saturated (e.g. dihydro) heterocycles count.
"""
from rdkit import Chem

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    
    A nucleoside should include:
      (1) A sugar ring: typically a five-membered furanose (non-aromatic) ring composed of exactly
          1 oxygen and 4 carbons.
      (2) A nucleobase: a heterocyclic ring that contains at least one nitrogen.
          Ideally it is aromatic, but we relax that to allow dihydro forms.
      (3) An N-glycosidic bond: a covalent bond connecting the anomeric carbon (typically in the sugar)
          to a nucleobase nitrogen. We require that the N is a member of a ring (of size 5 or 6) which 
          is different from the sugar ring.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule fits as a nucleoside, False otherwise.
        str: Human-readable reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ringInfo = mol.GetRingInfo()
    all_rings = ringInfo.AtomRings()
    
    # 1. Locate the sugar ring candidate: a non-aromatic 5-membered ring with exactly 1 oxygen and 4 carbons.
    sugar_ring = None
    for ring in all_rings:
        if len(ring) != 5:
            continue
        # Get the atoms in the ring
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # For sugars we expect none of the atoms to be flagged aromatic.
        if any(atom.GetIsAromatic() for atom in atoms):
            continue  # skip if any aromatic flag seen
        num_oxygen = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
        num_carbon = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        if num_oxygen == 1 and num_carbon == 4:
            sugar_ring = set(ring)
            break
    if sugar_ring is None:
        return False, "No ribose or deoxyribose sugar ring found (expected non-aromatic 5-membered ring with 1 oxygen and 4 carbons)"
    
    # 2. Identify candidate glycosidic bonds
    # Look for a bond where one end (in sugar_ring) is connected to a nitrogen (atomic number 7)
    # that is not within the sugar. In addition, we require that the nitrogen is part of a heterocyclic ring.
    glyco_bond_found = False
    nucleobase_ring_found = False
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        idx1 = a1.GetIdx()
        idx2 = a2.GetIdx()
        # Check if exactly one atom is in the sugar ring
        if (idx1 in sugar_ring and idx2 not in sugar_ring) or (idx2 in sugar_ring and idx1 not in sugar_ring):
            # Identify which is the sugar atom and which is the candidate nucleobase atom.
            sugar_atom = a1 if idx1 in sugar_ring else a2
            cand_atom = a2 if idx1 in sugar_ring else a1
            if cand_atom.GetAtomicNum() != 7:
                continue  # not a nitrogen, skip
            # Optionally, we can require that the sugar_atom is a carbon (the anomeric carbon should be C).
            if sugar_atom.GetAtomicNum() != 6:
                continue
            # Now check that the candidate nitrogen is part of a ring (of size 5 or 6) that is not the sugar ring.
            cand_in_rings = []
            for ring in all_rings:
                if cand_atom.GetIdx() in ring and not set(ring).issubset(sugar_ring):
                    cand_in_rings.append(ring)
            if not cand_in_rings:
                continue  # candidate nitrogen is not in a separate ring; skip it.
            # For at least one candidate ring, check that the ring is of appropriate size.
            valid_base = False
            for ring in cand_in_rings:
                if len(ring) in (5,6):
                    # As a minimal requirement, check that the ring has at least one other nitrogen.
                    num_nitrogens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
                    if num_nitrogens >= 1:  # allow even a single N, since the candidate itself is one.
                        valid_base = True
                        break
            if not valid_base:
                continue
            # If we passed all checks, we mark the glycosidic bond as found.
            glyco_bond_found = True
            nucleobase_ring_found = True
            break

    if not glyco_bond_found:
        return False, "No N-glycosidic bond between sugar and nucleobase detected"
    
    # 3. (Optional) Additional sanity check on molecular connectivity.
    # We could later add further tests (molecular weight, additional substituents), but here we assume that
    # the core criteria suffice.
    
    return True, "Molecule is a nucleoside: contains a sugar ring and a nucleobase (with N-heterocycle) linked via an N-glycosidic bond"

# Example usage:
if __name__ == "__main__":
    # A set of test SMILES (one example provided).
    test_smiles = "O1[C@@H](N2CN=CC=C2)[C@H](O)[C@H](O)[C@H]1CO"  # 1-(beta-D-Ribofuranosyl)-1,2-dihydropyrimidine
    result, reason = is_nucleoside(test_smiles)
    print(result, ":", reason)