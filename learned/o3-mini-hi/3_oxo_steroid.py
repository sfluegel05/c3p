"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: 3-oxo steroid
Definition (heuristic): A 3-oxo steroid is a molecule with a fused steroid nucleus
that resembles a cyclopentanoperhydrophenanthrene system – that is, four fused rings 
(typically, three six‐membered rings and one five‐membered ring or lactone) with 
roughly 15–21 carbon atoms – and that bears at least one ketone (C=O) group on one 
of the six‐membered rings.
Improvements over the previous attempt:
  - Require that the fused candidate steroid nucleus contains exactly 4 rings.
  - Accept a five–membered ring if it is either entirely carbons or is "lactone-like"
    (exactly 4 carbons and 1 oxygen).
  - Enforce a narrow range for the number of carbon atoms (15-21) in the nucleus.
  - Look only for a ketone on a six–membered ring within that nucleus.
"""
from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    
    A 3-oxo steroid is defined (heuristically) as a molecule that contains a fused
    steroid nucleus – the cyclopentanoperhydrophenanthrene system (4 fused rings with 
    about 15-21 carbon atoms, including at least one 5-membered ring or lactone) – 
    and that has at least one ketone (C=O) substituent on a six-membered ring of that nucleus.

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as a 3-oxo steroid, False otherwise.
        str : Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all rings (as tuples of atom indices) from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found; not a steroid nucleus"
    
    # Consider only rings of size 5 or 6.
    candidate_rings = [ring for ring in ring_info if len(ring) in (5, 6)]
    if len(candidate_rings) < 4:
        return False, "Fewer than 4 candidate rings (size 5 or 6), not a steroid nucleus"
    
    # Build a graph: each node is one candidate ring; rings are “fused” if they share at least 2 atoms.
    n = len(candidate_rings)
    ring_graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(set(candidate_rings[i]).intersection(candidate_rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
                
    # Find connected components (fused sets of rings)
    visited = set()
    components = []
    for i in range(n):
        if i not in visited:
            comp = set()
            stack = [i]
            while stack:
                cur = stack.pop()
                if cur not in comp:
                    comp.add(cur)
                    visited.add(cur)
                    stack.extend(ring_graph[cur] - comp)
            components.append(comp)
    
    # Choose the largest connected component as candidate steroid nucleus.
    largest_comp = max(components, key=lambda c: len(c))
    
    # Require that the fused nucleus consists of exactly 4 rings.
    if len(largest_comp) != 4:
        return False, "Fused ring system does not consist of exactly 4 rings typical of a steroid nucleus"
    
    # Collect all atom indices in this fused system.
    nucleus_atoms = set()
    for idx in largest_comp:
        nucleus_atoms.update(candidate_rings[idx])
    
    # Count carbon atoms in the nucleus.
    nucleus_carbon_count = sum(1 for atom in mol.GetAtoms() 
                                if atom.GetAtomicNum() == 6 and atom.GetIdx() in nucleus_atoms)
    if not (15 <= nucleus_carbon_count <= 21):
        return False, f"Fused ring nucleus has {nucleus_carbon_count} carbons; expected between 15 and 21"
    
    # Count the number of six-membered rings and check for an acceptable 5-membered ring.
    six_count = 0
    found_five = False  # will flag if a five-membered ring or lactone is found
    for idx in largest_comp:
        ring = candidate_rings[idx]
        if len(ring) == 6:
            six_count += 1
        elif len(ring) == 5:
            # Count element types in the ring.
            c_atoms = sum(1 for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 6)
            o_atoms = sum(1 for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 8)
            # Accept either a normal cyclopentane (5 carbons) or a lactone-like ring (4 carbons, 1 oxygen).
            if c_atoms == 5 or (c_atoms == 4 and o_atoms == 1):
                found_five = True
    if not found_five:
        return False, "Fused steroid nucleus does not contain a five-membered ring or lactone typical of steroids"
    
    # For typical steroids we expect 3 six-membered rings.
    required_six = 3 if nucleus_carbon_count >= 17 else 2
    if six_count < required_six:
        return False, f"Nucleus has only {six_count} six-membered ring(s); expected at least {required_six}"
    
    # Next, ensure that one of the six-membered rings in the nucleus carries a ketone (C=O)
    # Look for a carbon atom in one of these rings that is double-bonded to an oxygen.
    for idx in largest_comp:
        ring = candidate_rings[idx]
        if len(ring) == 6:
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() == 6:  # carbon
                    # Check each bond on this carbon.
                    for bond in atom.GetBonds():
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            nbr = bond.GetOtherAtom(atom)
                            if nbr.GetAtomicNum() == 8:
                                return True, ("Found ketone group on a six-membered ring of the steroid nucleus "
                                              "(likely at position 3)")
    
    return False, "No ketone group on any six-membered ring of the fused steroid nucleus was found"


# Sample test cases for demonstration:
if __name__ == "__main__":
    test_smiles = [
      # True examples (some from the provided list)
      "[H][C@@]12CCC3=CC(=O)C=C[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])C[C@@H](C(=O)OC)[C@]1(O)C(=O)CO",  # methyl prednisolone-16alpha-carboxylate
      "C[C@]12CC[C@H]3[C@@H](CC=C4CC(=O)CC[C@]34C)[C@@H]1CCC2=O",  # androst-5-ene-3,17-dione
      "[H][C@@]12CCC3=CC(=O)C=C[C@]3(C)[C@@]1([H])CC[C@]1(C)OC(=O)CC[C@@]21[H]",  # testolactone (lactone D-ring)
      # A false positive candidate (should return False)
      "C[C@H](CCC(=O)NCCCC(O)=O)[C@H]1CC[C@H]2[C@H]3[C@H](CC[C@]12C)[C@@]4(C)CC[C@@H](O)C[C@H]4CC3=O"
    ]
    for sm in test_smiles:
        result, reason = is_3_oxo_steroid(sm)
        print(sm.split()[0][:70] + " ...", "->", result, "|", reason)