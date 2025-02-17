"""
Classifies: CHEBI:73011 germacranolide
"""
"""
Classifies: Germacranolide – a sesquiterpene lactone based on the germacrane skeleton.
Improved criteria:
  (1) The SMILES must be parsed properly.
  (2) The molecule must contain a lactone ring (cyclic ester) detected with a SMARTS pattern.
  (3) The molecule must contain at least one macrocycle (or fused “ring system”) that:
         • has a total size typically between 10 and 12 atoms,
         • is predominantly carbon (>=80% C),
         • has at least one non‐aromatic sp2 carbon (indicating unsaturation),
         • and importantly, is fused with the lactone ring (i.e. the lactone atoms are part of this macrocycle).
  (4) The overall molecular weight should be in the typical range for sesquiterpene lactones (150–600 Da).
  (5) The total number of carbons should be near that of a sesquiterpene (15–22 carbons).
Note: This heuristic approach remains imperfect.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines whether the molecule given by the SMILES string is a germacranolide 
    (a sesquiterpene lactone based on a germacrane skeleton).

    Args:
        smiles (str): SMILES representation of the molecule.

    Returns:
        bool: True if the molecule is classified as a germacranolide, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Criterion (1): Look for the lactone ring (cyclic ester).
    # The SMARTS "[#6;R](=O)[OX2;R]" looks for a ring carbonyl carbon with an oxygen (both being in a ring).
    lactone_smarts = "[#6;R](=O)[OX2;R]"
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    if not lactone_matches:
        return False, "No lactone ring (cyclic ester) found"

    # Get all rings present in the molecule. Each ring is represented as a set of atomic indices.
    ssr_rings = [set(ring) for ring in mol.GetRingInfo().AtomRings()]
    if not ssr_rings:
        return False, "No rings detected in the molecule"

    # Merge overlapping rings to capture any fused ring systems.
    merged_rings = []
    while ssr_rings:
        current = ssr_rings.pop()
        merged = True
        while merged:
            merged = False
            remove_list = []
            for other in ssr_rings:
                if current.intersection(other):
                    current = current.union(other)
                    remove_list.append(other)
                    merged = True
            for rem in remove_list:
                ssr_rings.remove(rem)
        merged_rings.append(current)

    # Search for candidate macrocycle that meets the criteria.
    candidate_found = False
    candidate_reason = ""
    for ring_set in merged_rings:
        ring_size = len(ring_set)
        # We target macrocycles roughly of size 10 to 12 atoms.
        if not (10 <= ring_size <= 12):
            continue

        # Check that the ring is predominantly carbon (>=80%).
        carbon_count = sum(1 for idx in ring_set if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if carbon_count / ring_size < 0.8:
            continue

        # Check for unsaturation: at least one non-aromatic sp2 carbon within the ring.
        unsaturation = False
        for idx in ring_set:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2 and not atom.GetIsAromatic():
                unsaturation = True
                break
        if not unsaturation:
            continue

        # Additional key improvement:
        # The lactone ring should be fused with this candidate macrocyclic system.
        # We require that at least one lactone match (which is a tuple of atom indices) is entirely contained in the candidate ring.
        lactone_fused = False
        for match in lactone_matches:
            # We expect the lactone SMARTS to match at least two atoms in the ring (the carbonyl carbon and the ring oxygen).
            if set(match).issubset(ring_set):
                lactone_fused = True
                break
        if not lactone_fused:
            continue

        # If we reach here, candidate ring is appropriate.
        candidate_found = True
        candidate_reason = ("Contains a lactone ring fused with a macrocyclic (10–12 atoms, >=80% C) system "
                            "featuring unsaturation typical of a germacrane skeleton")
        break

    if not candidate_found:
        return False, ("No fused macrocyclic (10–12-membered) ring system with sufficient carbon content, unsaturation, "
                       "and association with the lactone ring was found")

    # Criterion (4): Check that overall molecular weight is in the range typical for sesquiterpene lactones.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150 or mol_wt > 600:
        return False, "Molecular weight not in range for typical sesquiterpene lactones (150–600 Da)"

    # Criterion (5): Check that total number of carbons is in the expected range for a sesquiterpene (15–22).
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 15 or total_carbons > 22:
        return False, f"Total number of carbons ({total_carbons}) is not in the expected range (15–22) for sesquiterpene lactones"

    return True, candidate_reason

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Provincialin is a known germacranolide example.
    test_smiles = r"C\C=C(/CO)C(=O)OC\C(=C/CO)C(=O)O[C@@H]1C\C(C)=C\C[C@H](OC(C)=O)\C(C)=C/[C@H]2OC(=O)C(=C)[C@H]12"
    result, reason = is_germacranolide(test_smiles)
    print(f"Input SMILES: {test_smiles}")
    print(f"Result: {result}\nReason: {reason}")