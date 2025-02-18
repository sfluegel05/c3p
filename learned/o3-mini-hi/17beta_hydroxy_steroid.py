"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
"""
Classifies: 17β-hydroxy steroid
Definition: A 17β-hydroxy steroid is defined as a steroid containing a fused carbocyclic nucleus with 4 rings – three 6‐membered (A, B, C rings) and one 5‐membered (D ring) – where one carbon in the five‐membered D–ring carries an –OH group at C17 with explicitly assigned (likely beta) stereochemistry.
The algorithm:
  1. Parse the SMILES.
  2. Identify candidate rings (5- or 6-membered mostly made of carbons).
  3. Build a “fusion graph” (rings fused if they share at least 2 atoms) and collect connected ring sets.
  4. Find a fused set with exactly 4 rings (3 six‐membered and 1 five‐membered), i.e. a steroid nucleus.
  5. From the nucleus, isolate the five–membered ring (assumed D–ring).
  6. For each carbon atom in the five–membered ring that is sp³, examine its neighbors.
     To be considered a candidate for the 17β –OH position, the carbon must have a single-bonded oxygen that carries at least one hydrogen, and aside from that oxygen (and hydrogens) every other heavy–atom neighbor must belong to the steroid nucleus. This is meant to ensure that the –OH is “integral” to the steroid core (and not part of a sugar or side–chain).
  7. If the candidate carbon has an explicitly set chiral tag then we return True with a note that the beta configuration is likely, otherwise the configuration is uncertain.
If any step fails, we return False with a reason.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines whether a molecule is a 17β-hydroxy steroid based on its SMILES string.
    
    The procedure is:
      (a) Parse the SMILES string.
      (b) Identify rings of size 5 or 6 that are mostly made of carbon:
          for 5-membered: require at least 4 carbons;
          for 6-membered: require at least 5 carbons.
      (c) Build a fusion graph among these candidate rings (two rings are fused if they share 2 or more atoms).
      (d) Find a fused component (connected set of rings) with exactly 4 rings: 3 six-membered and 1 five-membered.
      (e) From this nucleus, isolate the five-membered ring (the D-ring). Then, for each carbon (atom) from
          that ring (only sp3 carbons are considered) check for an attached -OH group.
          In addition, verify that aside from the oxygen the carbon is bonded only to atoms in the steroid nucleus.
      (f) If a candidate is found, return True along with an explanation indicating whether the –OH-bearing carbon
          has an explicit chiral tag (likely beta) or not.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): Tuple with True if the molecule appears to be a 17β-hydroxy steroid; otherwise False and an explanation.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings found in molecule"
    
    # Filter candidate rings: only those of size 5 or 6 which are largely carbocyclic.
    candidate_rings = []
    for ring in atom_rings:
        ring_size = len(ring)
        if ring_size not in (5, 6):
            continue
        n_carbons = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if ring_size == 5 and n_carbons < 4:
            continue
        if ring_size == 6 and n_carbons < 5:
            continue
        candidate_rings.append(set(ring))
    
    if not candidate_rings:
        return False, "No 5- or 6-membered (mostly) carbocyclic rings found (steroid nucleus expected)"
    
    # Build a fusion graph among candidate rings: rings are fused if they share 2 or more atoms.
    n_rings = len(candidate_rings)
    fusion = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if len(candidate_rings[i].intersection(candidate_rings[j])) >= 2:
                fusion[i].add(j)
                fusion[j].add(i)
    
    # Retrieve connected components in the fusion graph.
    visited = set()
    fused_components = []
    for i in range(n_rings):
        if i in visited:
            continue
        stack = [i]
        comp = set()
        while stack:
            cur = stack.pop()
            if cur in comp:
                continue
            comp.add(cur)
            visited.add(cur)
            for neigh in fusion[cur]:
                if neigh not in comp:
                    stack.append(neigh)
        fused_components.append(comp)
    
    # Look for a fused component that has exactly 4 rings:
    # count 3 six-membered rings and 1 five-membered ring.
    steroid_component = None
    for comp in fused_components:
        count6 = sum(1 for idx in comp if len(candidate_rings[idx]) == 6)
        count5 = sum(1 for idx in comp if len(candidate_rings[idx]) == 5)
        if count6 == 3 and count5 == 1:
            steroid_component = comp
            break
    if steroid_component is None:
        return False, "Steroid nucleus (3 six-membered and 1 five-membered fused rings) not found"
    
    # Collect all atom indices from the steroid nucleus.
    nucleus_atoms = set()
    five_membered_rings = []
    for idx in steroid_component:
        ring = candidate_rings[idx]
        nucleus_atoms.update(ring)
        if len(ring) == 5:
            five_membered_rings.append(ring)
    if not five_membered_rings:
        return False, "No 5-membered ring found in steroid nucleus (expected D-ring)"
    
    # Now examine each atom in the five-membered (D) ring.
    # Check for an -OH group attached to a sp3 carbon that is integral to the steroid nucleus.
    # "Integral" means that aside from the oxygen (and implicit hydrogens) every other heavy neighbor must be in the nucleus.
    for ring in five_membered_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or atom.GetHybridization() != rdchem.HybridizationType.SP3:
                continue
            # Check that aside from the oxygen (for -OH) all heavy neighbors belong to the nucleus.
            outside_neighbors = False
            candidate_O = None
            for nbr in atom.GetNeighbors():
                # Skip hydrogens (they might be implicit)
                if nbr.GetAtomicNum() == 1:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Look for an oxygen candidate that is single-bonded and carries at least one attached hydrogen.
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond is None or bond.GetBondTypeAsDouble() > 1.0:
                        continue
                    if nbr.GetTotalNumHs() < 1:
                        continue
                    candidate_O = nbr.GetIdx()
                else:
                    # If a heavy neighbor is not in the steroid nucleus, then it might be a substituent.
                    if nbr.GetIdx() not in nucleus_atoms:
                        outside_neighbors = True
                        break
            if outside_neighbors or candidate_O is None:
                continue
            # A candidate -OH on an sp3 carbon in the 5-membered ring (D-ring) is found.
            if atom.GetChiralTag() != rdchem.ChiralType.CHI_UNSPECIFIED:
                reason = ("Steroid nucleus (3 six-membered and 1 five-membered fused rings) found and an -OH "
                          "group is attached on a chiral (likely beta) carbon in the five-membered D-ring "
                          "consistent with a 17β–OH substituent.")
                return True, reason
            else:
                reason = ("Steroid nucleus (3 six-membered and 1 five-membered fused rings) found and an -OH "
                          "group is attached on a carbon in the five-membered D-ring; however, stereochemistry "
                          "is undefined (thus beta configuration is uncertain).")
                return True, reason
    return False, ("No -OH group attached on a steroid-core carbon of the five-membered ring was found "
                   "(expected for 17β–OH substitution)")

# Example usage: (Run as script to test with one positive example)
if __name__ == "__main__":
    # Example positive case: kidjoranin-3-O-beta-digitoxopyranoside
    test_smiles = "[H][C@@]1(C[C@H](O)[C@H](O)[C@@H](C)O1)O[C@H]1CC[C@@]2(C)C(C1)=CC[C@]1(O)[C@]2([H])C[C@@H](OC(=O)\\C=C\\c2ccccc2)[C@]2(C)[C@@](O)(CC[C@]12O)C(C)=O"
    result, explanation = is_17beta_hydroxy_steroid(test_smiles)
    print(result, explanation)