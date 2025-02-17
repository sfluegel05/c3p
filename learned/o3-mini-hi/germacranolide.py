"""
Classifies: CHEBI:73011 germacranolide
"""
"""
Classifies: Germacranolide – a sesquiterpene lactone based on the germacrane skeleton.

Improved criteria:
  (1) The SMILES must be parsed properly.
  (2) The molecule must contain a lactone ring (cyclic ester, detected with a SMARTS pattern).
  (3) The molecule must contain at least one macrocycle that is either a single ring
      or a fused “ring system” (a union of intersecting rings) with:
         • a total size between 10 and 12 atoms,
         • at least 80% of the atoms are carbons,
         • and at least one non‐aromatic sp2 carbon (sign of unsaturation) in that macrocycle.
  (4) The overall molecular weight should be in the typical range for sesquiterpene lactones (150-600 Da).
  (5) The total number of carbons is expected to be near that of a sesquiterpene (e.g. 15–22).

Note: This heuristic approach is still not perfect, but has been modified based on the
    results from the previous attempt.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines whether the molecule given by the SMILES string is a germacranolide (a sesquiterpene lactone
    based on a germacrane skeleton).

    Args:
        smiles (str): SMILES representation of the molecule.

    Returns:
        bool: True if the molecule is classified as a germacranolide, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Criterion (1): Check for lactone ring (cyclic ester).
    # This SMARTS looks for a carbonyl carbon in a ring bound to an oxygen (also in a ring).
    lactone_smarts = "[#6;R](=O)[OX2;R]"
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring (cyclic ester) found"

    # Get all SSSR rings as sets of atom indices.
    ssr_rings = [set(ring) for ring in mol.GetRingInfo().AtomRings()]
    if not ssr_rings:
        return False, "No rings detected in molecule"

    # Merge rings that overlap to capture fused ring systems.
    # We build connected components where rings sharing at least one atom are merged.
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

    # Now check each candidate macrocycle.
    candidate_found = False
    candidate_reason = ""
    for ring_set in merged_rings:
        ring_size = len(ring_set)
        # We want macrocycles of size 10 to 12.
        if not (10 <= ring_size <= 12):
            continue

        # Count carbons in the candidate ring.
        carbon_count = sum(1 for idx in ring_set if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if carbon_count / ring_size < 0.8:
            continue  # not predominantly carbon

        # Check for unsaturation within the fused system.
        # We require that within the atoms of this candidate ring, at least one carbon is sp2
        # (and not aromatic) suggesting a double (or equivalent) bond.
        unsaturation = False
        for idx in ring_set:
            atom = mol.GetAtomWithIdx(idx)
            # Only consider carbon atoms.
            if atom.GetAtomicNum() != 6:
                continue
            # Check if sp2 hybridized and non-aromatic.
            if atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2 and not atom.GetIsAromatic():
                unsaturation = True
                break

        if not unsaturation:
            continue

        # If one candidate ring system fulfills conditions, accept.
        candidate_found = True
        candidate_reason = "Contains a lactone ring and a fused macrocyclic (10–12 atoms, >=80% C) system with unsaturation typical of a germacrane skeleton"
        break

    if not candidate_found:
        return False, "No 10–12-membered ring (or fused ring system) with sufficient carbon content and unsaturation found"

    # Criterion (4): Check overall molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150 or mol_wt > 600:
        return False, "Molecular weight not in range for typical sesquiterpene lactones (150–600 Da)"

    # Criterion (5): Check total number of carbons.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 15 or total_carbons > 22:
        return False, f"Total number of carbons ({total_carbons}) is not in the expected range (15–22) for sesquiterpene lactones"

    return True, candidate_reason

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Example: Provincialin (should be classified as germacranolide)
    test_smiles = r"C\C=C(/CO)C(=O)OC\C(=C/CO)C(=O)O[C@@H]1C\C(C)=C\C[C@H](OC(C)=O)\C(C)=C/[C@H]2OC(=O)C(=C)[C@H]12"
    result, reason = is_germacranolide(test_smiles)
    print(f"Input SMILES: {test_smiles}")
    print(f"Result: {result}\nReason: {reason}")