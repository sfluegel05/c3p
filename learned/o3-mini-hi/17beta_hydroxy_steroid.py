"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
"""
Classifies: 17β-hydroxy steroid
Definition: A 17β-hydroxy steroid is defined as a steroid (with a fused nucleus of 4 rings – three six‐membered and one five‐membered, all carbocyclic)
            that carries a –OH substituent at the 17–position on the five–membered (D) ring with a defined (beta-like) stereochemistry.
            
Improved approach:
  1. Parse the SMILES.
  2. From the molecule’s ring information, collect only rings that are 5– or 6–membered and made entirely of carbons.
  3. Build a connectivity graph among these rings (two rings are “fused” if they share two or more atoms).
  4. Look for a fused component that has exactly 4 rings (3 of size 6 and 1 of size 5). This is our steroid nucleus.
  5. Then, from that nucleus, select the unique five–membered ring and search its carbon atoms for a substituent –OH.
  6. For any such candidate the carbon must be sp3 and carry an explicit chiral tag – if so, assume it is a 17β–OH.
  7. If an –OH arm is found on the five–membered ring but its chirality is undefined, return a result that shows the configuration is uncertain.
  
Note: This heuristic may miss some steroids and may still yield some mis‐classifications.
"""

from rdkit import Chem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines whether a molecule is a 17β-hydroxy steroid based on its SMILES string.
    
    The procedure is:
      (a) Parse SMILES; if none then return False.
      (b) Identify rings that are either 5– or 6–membered and are completely carbocyclic.
      (c) Build a fusion graph among these rings (rings are fused if they share 2 or more atoms).
      (d) Look for a fused system that contains exactly 4 rings: three 6–membered and one 5–membered.
      (e) In that fused system, isolate the 5–membered ring (steroid D-ring) and search its member carbons
          for an –OH substituent. For each, require that the carbon is sp3 and has an explicitly assigned chirality.
      (f) If such a –OH group is found, return True along with an explanation.
          If the –OH group is found on the 5–membered ring but lacks chiral tag, return True but explain that the beta configuration is uncertain.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): A tuple where the first element is True if the molecule appears to be a 17β–hydroxy steroid,
                   and False otherwise; the second element is a string that explains the reasoning.
                   If SMILES parsing fails, returns (False, "Invalid SMILES string").
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings found in molecule"
    
    # Only consider rings that are 5 or 6 members AND with all atoms being carbon.
    candidate_rings = []
    for ring in atom_rings:
        if len(ring) not in (5, 6):
            continue
        if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            candidate_rings.append(set(ring))
    if not candidate_rings:
        return False, "No 5- or 6-membered pure-carbon rings found (steroid nucleus expected to be carbocyclic)"
    
    # Build a graph among these candidate rings: two rings are fused if they share 2 or more atoms.
    n_rings = len(candidate_rings)
    fusion = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if len(candidate_rings[i].intersection(candidate_rings[j])) >= 2:
                fusion[i].add(j)
                fusion[j].add(i)
    
    # Find connected components in the fusion graph.
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
    
    # Look for a fused component that has exactly 4 rings:  three six-membered and one five-membered.
    steroid_component = None
    for comp in fused_components:
        count6 = 0
        count5 = 0
        for idx in comp:
            ring_size = len(candidate_rings[idx])
            if ring_size == 6:
                count6 += 1
            elif ring_size == 5:
                count5 += 1
        if count6 == 3 and count5 == 1:
            steroid_component = comp
            break
    if steroid_component is None:
        return False, "Steroid nucleus (exactly 3 six-membered and 1 five-membered fused carbocyclic rings) not found"
    
    # Get all the atom indices in the steroid nucleus.
    nucleus_atoms = set()
    five_membered_rings = []  # there should be one but we make a list in case
    for idx in steroid_component:
        ring = candidate_rings[idx]
        nucleus_atoms.update(ring)
        if len(ring) == 5:
            five_membered_rings.append(ring)
    
    if not five_membered_rings:
        return False, "No 5-membered ring found in steroid nucleus (expected for the D-ring)"
    
    # Now, examine atoms in the five-membered ring to see if one bears an -OH group.
    # Instead of relying solely on a SMARTS, we loop through the atoms in the five-membered ring
    # and check for a directly attached oxygen that also carries at least one hydrogen (implying -OH).
    for ring in five_membered_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Only consider sp3 carbons (to ensure aliphatic, non-carbonyl context)
            if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
            # Look at neighbors
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    # Check that this oxygen is –OH. In RDKit hydrogens may be implicit.
                    # We require that it has at least one attached hydrogen either implicitly or explicitly.
                    numHs = nbr.GetTotalNumHs()
                    if numHs < 1:
                        continue
                    # We also want to ensure that the oxygen is not part of a carbonyl. A carbonyl oxygen will be double-bonded.
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond is not None and bond.GetBondTypeAsDouble() > 1.0:
                        continue
                    # Candidate found. Now check chirality on the carbon.
                    if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                        reason = ("Steroid nucleus (3 six-membered, 1 five-membered, all carbocyclic) found and an -OH "
                                  "group is attached on a chiral (likely beta) carbon in the five-membered D ring "
                                  "consistent with a 17β–OH substituent.")
                        return True, reason
                    else:
                        reason = ("Steroid nucleus (3 six-membered, 1 five-membered, all carbocyclic) found and an -OH "
                                  "group is attached in the five-membered D ring; however, the stereochemistry at that carbon "
                                  "is not defined (thus beta configuration is uncertain).")
                        return True, reason
    return False, "No -OH group found on any carbon of the 5-membered ring in the steroid nucleus (expected for 17β–OH)"


# Example usage:
if __name__ == "__main__":
    # Example positive case: (kidjoranin-3-O-beta-digitoxopyranoside contains a steroid nucleus with an -OH in the D ring)
    test_smiles = "[H][C@@]1(C[C@H](O)[C@H](O)[C@@H](C)O1)O[C@H]1CC[C@@]2(C)C(C1)=CC[C@]1(O)[C@]2([H])C[C@@H](OC(=O)\\C=C\\c2ccccc2)[C@]2(C)[C@@](O)(CC[C@]12O)C(C)=O"
    result, explanation = is_17beta_hydroxy_steroid(test_smiles)
    print(result, explanation)