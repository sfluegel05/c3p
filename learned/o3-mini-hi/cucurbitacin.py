"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: Cucurbitacin (tetracyclic triterpenoids from the cucurbitane family)
This improved classifier first selects the likely “aglycone” (by choosing the fragment with the most carbon atoms)
before assessing:
  – A fused ring system of at least 4 rings,
  – At least 27 carbon atoms (triterpenoid skeleton),
  – Presence of at least one carbonyl (C=O).
Note: Cucurbitacins are often decorated with sugars and extra oxygens. This heuristic‐based classifier
may miss or mis‐classify edge cases.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def get_aglycone(mol):
    """
    Given a (possibly glycosylated) molecule, break it into fragments
    and return the fragment with the highest number of carbon atoms.
    """
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) == 1:
        return mol
    best_frag = None
    max_carbons = 0
    for frag in frags:
        carbon_count = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
        if carbon_count > max_carbons:
            max_carbons = carbon_count
            best_frag = frag
    return best_frag if best_frag is not None else mol

def has_fused_tetracyclic_core(mol):
    """
    Checks whether the molecule has at least one connected set of 4 (or more) fused rings.
    Each ring is represented as a set of atom indices (from Mol.GetRingInfo().AtomRings()).
    Two rings are considered connected if they share at least one atom.
    Then we group rings into connected components. One component must have at least 4 rings.
    """
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) < 4:
        return False, f"Found only {len(atom_rings)} rings in the aglycone; expecting at least 4 fused rings."
    
    # Build a graph where nodes represent rings and an edge exists if two rings share any atom.
    ring_sets = [set(r) for r in atom_rings]
    n = len(ring_sets)
    graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if ring_sets[i].intersection(ring_sets[j]):
                graph[i].add(j)
                graph[j].add(i)
    # Now find connected components:
    seen = set()
    components = []
    for i in range(n):
        if i not in seen:
            stack = [i]
            comp = set()
            while stack:
                node = stack.pop()
                if node in comp: 
                    continue
                comp.add(node)
                for neigh in graph[node]:
                    if neigh not in comp:
                        stack.append(neigh)
            seen = seen.union(comp)
            components.append(comp)
    # Check if any connected component has at least 4 rings.
    for comp in components:
        if len(comp) >= 4:
            return True, ("Aglycone fragment shows a fused ring system "
                          f"with {len(comp)} rings, consistent with a tetracyclic core.")
    return False, "No connected set of at least 4 fused rings found in the aglycone fragment."

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids (often bearing sugar moieties) and typically
    display a fused tetracyclic core (~30 carbons) and one or more ketone groups.
    
    The procedure is:
      1. Parse the SMILES string.
      2. From the possibly multi-fragment molecule, obtain the aglycone (largest carbon count).
      3. Check that the aglycone has a fused ring system of at least 4 rings.
      4. Check that the aglycone has at least 27 carbon atoms.
      5. Check for the presence of one or more carbonyl groups (C=O) in the aglycone.
      6. Optionally, check that the molecular weight of the aglycone is in the expected range (>=400 Da).
    
    Returns:
        (bool, str): A tuple (True, reason) if classified as cucurbitacin; (False, reason) otherwise.
    """
    # Parse
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Get aglycone fragment (aiming to remove sugar moieties)
    aglycone = get_aglycone(mol)
    
    # Check fused ring system
    fused_ok, ring_reason = has_fused_tetracyclic_core(aglycone)
    if not fused_ok:
        return False, ring_reason

    # Count carbons in aglycone
    carbon_count = sum(1 for atom in aglycone.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 27:
        return False, f"Aglycone has only {carbon_count} carbon atoms; expected at least 27 for a triterpenoid scaffold."
    
    # Check for at least one carbonyl group (ketone)
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=O")
    carbonyl_matches = aglycone.GetSubstructMatches(carbonyl_pattern)
    if len(carbonyl_matches) < 1:
        return False, "No carbonyl group (C=O) detected in the aglycone; cucurbitacins typically feature ketone functionality."
    
    # Optional: Check molecular weight of the aglycone.
    mol_wt = rdMolDescriptors.CalcExactMolWt(aglycone)
    if mol_wt < 400:
        return False, f"Aglycone molecular weight is too low ({mol_wt:.1f} Da) compared with typical cucurbitacin aglycones."
    
    return True, ("Aglycone fragment exhibits a fused tetracyclic system with sufficient carbon count "
                  f"({carbon_count} C atoms) and at least one carbonyl group, and has a molecular weight of {mol_wt:.1f} Da. "
                  "These features are consistent with a cucurbitacin (tetracyclic triterpenoid) scaffold.")

# Example usage:
if __name__ == "__main__":
    # Example with cucurbitacin D (one of the provided examples)
    example_smiles = "[H][C@@]12C[C@H](O)C(=O)C(C)(C)C1=CC[C@@]1([H])[C@]3(C)C[C@@H](O)[C@]([H])([C@@](C)(O)C(=O)\\C=C\\C(C)(C)O)[C@@]3(C)CC(=O)[C@@]21C"
    result, reason = is_cucurbitacin(example_smiles)
    print("Cucurbitacin check:", result)
    print("Reason:", reason)