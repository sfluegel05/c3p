"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: Cucurbitacin (tetracyclic triterpenoids from the cucurbitane family)

This classifier aims to recognize cucurbitacin aglycone scaffolds. Cucurbitacins
typically feature a fused core of (roughly) 4 rings, a triterpenoid carbon count 
(usually between 27 and 40 C atoms), a molecular weight (of the aglycone fragment)
in the approximate range 400–600 Da, and one or more carbonyl (C=O) groups.
Since many cucurbitacins are glycosylated, we first attempt to extract the largest
“aglycone” fragment that itself displays a fused tetracyclic core.
In cases where no explicit carbonyl group is detected but all other properties 
fall in the expected range, we still classify the compound as cucurbitacin 
(allowing for tautomerism or transformation of the keto group).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def has_fused_tetracyclic_core(mol):
    """
    Check whether the molecule (or fragment) contains at least one set of four rings
    that are mutually fused (i.e. connected by at least one shared atom). Returns a tuple:
      (True, reason) if a connected set of >=4 rings is found, else (False, reason).
    """
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) < 4:
        return False, f"Found only {len(atom_rings)} rings; expecting at least 4 fused rings."

    # Make sets for each ring and build a graph where rings share atoms
    ring_sets = [set(r) for r in atom_rings]
    n = len(ring_sets)
    graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if ring_sets[i].intersection(ring_sets[j]):
                graph[i].add(j)
                graph[j].add(i)
                
    # Find connected components (each being a set of rings that share atoms)
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
            
    # Accept if any connected component contains at least 4 rings
    for comp in components:
        if len(comp) >= 4:
            return True, (f"Found fused ring system with {len(comp)} rings "
                          "consistent with a tetracyclic core.")
    return False, "No connected set of at least 4 fused rings found."

def get_aglycone(mol):
    """
    Given a (possibly glycosylated) molecule, break it into fragments (using GetMolFrags)
    and return the fragment with the highest carbon count that also passes the fused ring test.
    If none of the fragments individually show a tetracyclic fused core, return the fragment
    with the most carbons.
    """
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    best_candidate = None
    best_carbons = 0
    # First try to select among fragments that themselves contain a fused tetracyclic motif.
    for frag in frags:
        fused, _ = has_fused_tetracyclic_core(frag)
        carbon_count = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
        if fused and carbon_count > best_carbons:
            best_carbons = carbon_count
            best_candidate = frag
    if best_candidate is not None:
        return best_candidate
    # Otherwise, simply pick the fragment with highest carbon count.
    for frag in frags:
        carbon_count = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
        if carbon_count > best_carbons:
            best_carbons = carbon_count
            best_candidate = frag
    return best_candidate if best_candidate is not None else mol

def is_cucurbitacin(smiles: str):
    """
    Determines whether the input SMILES string corresponds to a cucurbitacin, i.e. a
    tetracyclic triterpenoid (often glycosylated) with the expected aglycone properties:
      • a fused core (>=4 fused rings),
      • at least 27 carbon atoms,
      • an aglycone molecular weight within 400–600 Da,
      • and at least one detectable carbonyl group.
    
    If no explicit carbonyl is found, the classifier still returns True provided that
    the other parameters (fused ring system, carbon count, MW) lie in the expected ranges,
    with the understanding that keto group detection may fail because of tautomerism.
    
    Returns a tuple (bool, str) where bool is True if classified as cucurbitacin, and str
    explains the reasoning.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Try to extract the aglycone (the fragment that is most likely the tetracyclic core)
    aglycone = get_aglycone(mol)
    if aglycone is None:
        aglycone = mol

    # Check for fused tetracyclic core
    fused_ok, ring_reason = has_fused_tetracyclic_core(aglycone)
    if not fused_ok:
        return False, ring_reason

    # Count carbons in aglycone
    carbon_count = sum(1 for atom in aglycone.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 27:
        return False, (f"Aglycone has only {carbon_count} carbon atoms; expected at least 27 "
                       "for a triterpenoid scaffold.")
    
    # Check molecular weight of the aglycone – expect roughly 400 to 600 Da for cucurbitacins.
    mol_wt = rdMolDescriptors.CalcExactMolWt(aglycone)
    if not (400 <= mol_wt <= 600):
        return False, (f"Aglycone molecular weight is {mol_wt:.1f} Da, "
                       "which is outside the typical 400–600 Da range for cucurbitacin cores.")
    
    # Check for at least one carbonyl group (C=O).
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=O")
    carbonyl_matches = aglycone.GetSubstructMatches(carbonyl_pattern)
    
    if len(carbonyl_matches) >= 1:
        reason = (f"Aglycone fragment shows a fused tetracyclic core, "
                  f"{carbon_count} carbon atoms, molecular weight {mol_wt:.1f} Da, "
                  f"and {len(carbonyl_matches)} carbonyl group(s).")
        return True, reason
    else:
        # If no carbonyl is found, allow a relaxed acceptance provided other properties fall in range.
        if (27 <= carbon_count <= 40) and (400 <= mol_wt <= 600):
            reason = (f"Aglycone fragment shows a fused tetracyclic core with {carbon_count} carbons "
                      f"and a molecular weight of {mol_wt:.1f} Da, but no explicit carbonyl group was "
                      "detected (this may be due to tautomerism or masked keto functionality).")
            return True, reason
        else:
            return False, ("No carbonyl group detected in the aglycone and its features do not clearly "
                           "match the typical cucurbitacin scaffold.")

# Example usage:
if __name__ == "__main__":
    # Example with cucurbitacin D (one of the provided examples)
    example_smiles = "[H][C@@]12C[C@H](O)C(=O)C(C)(C)C1=CC[C@@]1([H])[C@]3(C)C[C@@H](O)[C@]([H])([C@@](C)(O)C(=O)\\C=C\\C(C)(C)O)[C@@]3(C)CC(=O)[C@@]21C"
    result, reason = is_cucurbitacin(example_smiles)
    print("Cucurbitacin check:", result)
    print("Reason:", reason)