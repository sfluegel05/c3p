"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: Flavones
A flavone is defined as any flavonoid having a 2-aryl-1-benzopyran-4-one (2-phenylchromen-4-one) skeleton.
This function tests whether a given SMILES string contains that core and whether the rest of the molecule
is consistent with typical flavonoids (e.g. substitutions are relatively small or sugar‐like rather than forming
an entirely separate large scaffold).
"""

from rdkit import Chem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone (or a substituted flavone derivative) based on its SMILES string.
    The criteria are:
      1. The molecule must contain a 2-phenylchromen-4-one (flavone) core as detected by a SMARTS pattern.
      2. The matching fragment should have an expected number of heavy atoms (about 15)
         and show three fused/aromatic ring systems (the benzopyranone core plus the pendant phenyl).
      3. External substituents off the core must not be unusually large (unless they are sugar-like).
         Here we check that any group attached to the core is small (≤10 heavy atoms) or looks like a sugar (5–7 heavy atoms with ≥2 oxygens).
         
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the flavone core is found and the substituents are consistent with expected flavonoids.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the 2-phenylchromen-4-one core.
    # This pattern is intended to capture a benzopyran-4-one fused system with an attached phenyl.
    # (Note: Substituents on the rings will not be part of the match.)
    flavone_pattern = Chem.MolFromSmarts("c1ccc2c(c1)oc(=O)c(c2)-c3ccccc3")
    if flavone_pattern is None:
        return False, "Error in SMARTS for flavone core"
    
    # Find all substructure matches of the flavone core.
    matches = mol.GetSubstructMatches(flavone_pattern)
    if not matches:
        return False, "Does not contain a 2-phenylchromen-4-one (flavone) core"
    
    # Helper: retrieve the connected component (set of atoms) starting from a seed atom, but not including any atom in a given set.
    def get_component(start_idx, exclude_set):
        comp = set()
        stack = [start_idx]
        while stack:
            current = stack.pop()
            if current in comp:
                continue
            comp.add(current)
            for nbr in mol.GetAtomWithIdx(current).GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx not in exclude_set and nbr_idx not in comp:
                    stack.append(nbr_idx)
        return comp
    
    # For each match candidate, apply additional checks.
    for match in matches:
        # We expect the flavone core fragment (if unsubstituted) to have 15 heavy atoms
        # (flavone formula is C15H10O2; here we count non-hydrogen atoms).
        if len(match) != 15:
            # The match might have extra atoms if substitutions are forcibly included
            # so we require exactly 15 atoms to ensure we capture only the core.
            continue

        # Extract the fragment as a sub-molecule.
        try:
            submol = Chem.PathToSubmol(mol, list(match))
        except Exception:
            continue

        # Check ring count in the substructure.
        # A typical flavone includes three rings: two fused (forming the benzopyranone) plus the phenyl ring.
        rings = submol.GetRingInfo().AtomRings()
        # Count only rings with 5-7 atoms (most aromatic rings fall in this range).
        ring_count = sum(1 for r in rings if 5 <= len(r) <= 7)
        if ring_count != 3:
            # This candidate does not show the expected set of three rings.
            continue

        # Now check the attachments: for atoms in the core, check external substituents.
        core_set = set(match)
        external_issue = False
        for idx in core_set:
            atom = mol.GetAtomWithIdx(idx)
            # For every neighbor that is outside the core, get its connected component
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in core_set:
                    continue
                comp = get_component(nbr_idx, core_set)
                # Count heavy atoms in the external fragment (exclude hydrogens)
                heavy_count = sum(1 for i in comp if mol.GetAtomWithIdx(i).GetAtomicNum() > 1)
                # Count oxygens in that fragment.
                oxy_count = sum(1 for i in comp if mol.GetAtomWithIdx(i).GetAtomicNum() == 8)
                # Heuristic: allow small substituents (≤10 heavy atoms) or sugar-like fragments:
                # if fragment is 5-7 heavy atoms and has at least 2 oxygens, accept it.
                if heavy_count > 10:
                    external_issue = True
                    break
                if heavy_count >= 8 and not (5 <= heavy_count <= 7 and oxy_count >= 2):
                    # if moderate in size but not sugar-like, also flag it.
                    external_issue = True
                    break
            if external_issue:
                break

        if external_issue:
            # This candidate shows extra large substituents attached to the core,
            # and so its overall architecture is not dominated by the typical flavone motif.
            continue

        # If we get here, then the candidate has:
        #  - exactly 15 atoms in the match,
        #  - exactly 3 aromatic rings in the core,
        #  - and only relatively small/sugar-like substituents attached.
        return True, "Contains a complete flavone (2-phenylchromen-4-one) core with acceptable substituents."
    
    # If no candidate passed all checks, report false.
    # (Note that if a candidate matched the simple SMARTS but failed extra checks, it is not classified as a flavone.)
    return False, "Flavone core detected but substituents or ring system are not consistent with typical flavones."

# For demonstration purposes only:
if __name__ == "__main__":
    # You may test with one or more SMILES strings.
    # For example, 3-methoxyapigenin (a flavone) SMILES:
    test_smiles = "COc1c(oc2cc(O)cc(O)c2c1=O)-c1ccc(O)cc1"
    res, reason = is_flavones(test_smiles)
    print("Test molecule:", test_smiles)
    print("Is flavone?", res)
    print("Reason:", reason)