"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
"""
Classifies: 17β-hydroxy steroid
Definition: A 17β-hydroxy steroid is defined as a steroid containing a fused nucleus composed of 4 rings – three 6‐membered and one 5‐membered – that is mostly carbocyclic, and having an –OH substituent at the 17–position on the five–membered (D) ring with defined beta stereochemistry (if assigned).

Improved approach:
  1. Parse the SMILES string.
  2. Identify rings (from RDKit ring info) that are either 5- or 6-membered and are nearly all carbons;
     for 5-membered rings require at least 4 carbons, and for 6-membered rings at least 5 carbons.
  3. Build a “fusion graph” among these rings; rings are considered fused if they share at least 2 atoms.
  4. Search the connected components for one that has exactly 4 rings: 3 six-membered and 1 five-membered.
  5. In that candidate steroid nucleus, isolate the unique 5-membered ring (the D-ring).
  6. For each atom in that ring (taking only sp³ carbons), look for a neighbor oxygen that carries at least one hydrogen 
     (indicating an –OH group). Check whether the carbon atom bearing the –OH has an explicitly assigned chiral tag.
  7. Return True along with an explanation if a candidate is found; if the stereochemistry is missing, return True with
     a note that the beta configuration is uncertain.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines whether a molecule is a 17β-hydroxy steroid based on its SMILES string.
    
    The procedure is:
      (a) Parse the SMILES string.
      (b) Identify rings that are either 5- or 6-membered and are largely carbocyclic:
          for rings of size 5, allow if at least 4 atoms are carbon;
          for rings of size 6, allow if at least 5 atoms are carbon.
      (c) Build a fusion graph among these candidate rings (two rings are declared fused if they share 2 or more atoms).
      (d) Find a connected set of rings (a fused component) that contains exactly 4 rings: three six-membered and one five-membered.
      (e) In that component, identify the five-membered ring (the D-ring) and search its member carbons for an –OH substituent.
           For each candidate, the carbon must be sp3; then if it carries an explicit chiral tag assume (likely beta) configuration.
      (f) If found, return True along with an explanation. If the –OH substituent is found but the chiral tag is missing,
          return True with a note that the beta configuration is uncertain.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): Tuple where the first element is True if the molecule appears to be a 17β-hydroxy steroid,
                   and False otherwise. The second element explains the reasoning. If parsing fails, returns (False, "Invalid SMILES string").
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings found in molecule"
    
    # Filter candidate rings: only consider rings that have 5 or 6 members and are mostly carbons.
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
    
    # Build a fusion graph among candidate rings: two rings are fused if they share 2 or more atoms.
    n_rings = len(candidate_rings)
    fusion = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if len(candidate_rings[i].intersection(candidate_rings[j])) >= 2:
                fusion[i].add(j)
                fusion[j].add(i)
    
    # Retrieve connected components from the fusion graph
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
    
    # Look for a component that has exactly 4 rings: 3 six-membered and 1 five-membered.
    steroid_component = None
    for comp in fused_components:
        count6 = sum(1 for idx in comp if len(candidate_rings[idx]) == 6)
        count5 = sum(1 for idx in comp if len(candidate_rings[idx]) == 5)
        if count6 == 3 and count5 == 1:
            steroid_component = comp
            break
    if steroid_component is None:
        return False, "Steroid nucleus (3 six-membered and 1 five-membered fused rings) not found"
    
    # Collect all candidate atom indices in the steroid nucleus and record the five-membered rings.
    nucleus_atoms = set()
    five_membered_rings = []
    for idx in steroid_component:
        ring = candidate_rings[idx]
        nucleus_atoms.update(ring)
        if len(ring) == 5:
            five_membered_rings.append(ring)
    if not five_membered_rings:
        return False, "No 5-membered ring found in steroid nucleus (expected D-ring)"
    
    # Examine each atom in the five-membered ring (D-ring) for an attached -OH group.
    for ring in five_membered_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Only consider carbon atoms with sp3 hybridization (i.e. aliphatic)
            if atom.GetAtomicNum() != 6 or atom.GetHybridization() != rdchem.HybridizationType.SP3:
                continue
            # For each neighbor, if it is oxygen then check if it is an -OH group.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() != 8:
                    continue
                # Check if the oxygen has at least one hydrogen (even if implicit),
                # and ensure it is not involved in a double bond (e.g. as in a carbonyl).
                numHs = nbr.GetTotalNumHs()
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if numHs >= 1 and bond is not None and bond.GetBondTypeAsDouble() <= 1.0:
                    # Candidate -OH found. Now check chirality on the carbon.
                    if atom.GetChiralTag() != rdchem.ChiralType.CHI_UNSPECIFIED:
                        reason = ("Steroid nucleus (3 six-membered and 1 five-membered fused rings) found and an -OH "
                                  "group is attached on a chiral (likely beta) carbon in the five-membered D-ring "
                                  "consistent with a 17β–OH substituent.")
                        return True, reason
                    else:
                        reason = ("Steroid nucleus (3 six-membered and 1 five-membered fused rings) found and an -OH "
                                  "group is attached on a carbon in the five-membered D-ring; however, the stereochemistry "
                                  "at that carbon is undefined (thus beta configuration is uncertain).")
                        return True, reason
    return False, "No -OH group found on any carbon of the 5-membered ring in the steroid nucleus (expected for 17β–OH)"

# Example usage:
if __name__ == "__main__":
    # Example positive case: kidjoranin-3-O-beta-digitoxopyranoside
    test_smiles = "[H][C@@]1(C[C@H](O)[C@H](O)[C@@H](C)O1)O[C@H]1CC[C@@]2(C)C(C1)=CC[C@]1(O)[C@]2([H])C[C@@H](OC(=O)\\C=C\\c2ccccc2)[C@]2(C)[C@@](O)(CC[C@]12O)C(C)=O"
    result, explanation = is_17beta_hydroxy_steroid(test_smiles)
    print(result, explanation)