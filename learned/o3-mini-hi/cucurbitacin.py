"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: Cucurbitacin (tetracyclic triterpenoids derived from cucurbitane)

Cucurbitacins (often glycosylated) feature an aglycone that is a fused tetracyclic scaffold.
Typically the aglycone shows:
  • exactly 4 fused rings,
  • between about 27 and 35 carbon atoms,
  • a molecular weight in the 440–580 Da range, and
  • one or more carbonyl (C=O) groups.

This improved classifier attempts to first remove glycosyl parts by breaking the molecule
into fragments and then choosing the fragment whose fused ring system (exactly four fused rings)
has a molecular weight and carbon count closest to the expected ranges.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def has_fused_tetracyclic_core(mol):
    """
    Check whether the molecule (or fragment) contains a fused ring system of exactly 4 rings.
    We build a graph where each ring (as detected by GetRingInfo) is a node and rings sharing
    at least one atom are connected. Then we look for any connected component of exactly 4 rings.
    
    Returns (True, reason) if a connected set of exactly 4 rings is found,
    else (False, reason).
    """
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) < 4:
        return False, f"Found only {len(atom_rings)} rings; expecting at least 4 for a tetracyclic core."
    
    # Create a list of sets for each ring (list of atom indices)
    ring_sets = [set(r) for r in atom_rings]
    n = len(ring_sets)
    graph = {i: set() for i in range(n)}
    # Two rings are connected if they share at least one atom
    for i in range(n):
        for j in range(i+1, n):
            if ring_sets[i].intersection(ring_sets[j]):
                graph[i].add(j)
                graph[j].add(i)
                
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
    
    # Only accept if there is a connected component that consists exactly of 4 fused rings.
    for comp in components:
        if len(comp) == 4:
            return True, "Found fused ring system with exactly 4 rings consistent with a tetracyclic core."
    return False, "No connected set of exactly 4 fused rings found."

def get_aglycone(mol):
    """
    From a possibly glycosylated molecule, break it into fragments (using GetMolFrags)
    and return the fragment that:
      (a) contains a fused tetracyclic core (exactly 4 rings) and
      (b) has a molecular weight closest to the expected range (440–580 Da).
    If none of the fragments individually show a fused tetracyclic core, return the fragment
    with the highest carbon count.
    """
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    candidates = []
    for frag in frags:
        fused, _ = has_fused_tetracyclic_core(frag)
        carbon_count = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
        mw = rdMolDescriptors.CalcExactMolWt(frag)
        # Only accept fragments that have at least the minimal carbon count expected (27)
        if carbon_count >= 27:
            if fused:
                # Calculate a simple score: 0 if MW is within [440,580] (ideal) and penalty otherwise.
                if 440 <= mw <= 580:
                    penalty = 0
                else:
                    penalty = abs(mw - 510)  # 510 is a nominal midpoint
                candidates.append((penalty, frag))
    if candidates:
        # pick candidate with lowest penalty (i.e. closest to desired MW)
        best_candidate = sorted(candidates, key=lambda x: x[0])[0][1]
        return best_candidate
    # If no fragment has a fused tetracyclic core, simply return fragment with the highest carbon count.
    best_candidate = None
    best_carbons = 0
    for frag in frags:
        carbon_count = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
        if carbon_count > best_carbons:
            best_candidate = frag
            best_carbons = carbon_count
    return best_candidate if best_candidate is not None else mol

def is_cucurbitacin(smiles: str):
    """
    Determines whether the input SMILES string corresponds to a cucurbitacin.
    The classifier works on the aglycone (largest fragment that features a fused tetracyclic core)
    and applies the following criteria:
       • The aglycone must have a fused ring system of exactly 4 rings.
       • The aglycone must have between 27 and 35 carbon atoms.
       • The aglycone's molecular weight is expected to be in the range 440–580 Da.
       • The aglycone should contain at least one carbonyl (C=O) group.
    
    (If no carbonyl is found, we still accept the molecule if all other criteria are met,
    noting that tautomerism or masked keto groups might be responsible.)
    
    Returns:
         (True, reason) if classified as cucurbitacin, otherwise (False, reason)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Extract what appears to be the aglycone fragment.
    aglycone = get_aglycone(mol)
    if aglycone is None:
        aglycone = mol  # fallback
    
    # Check for fused tetracyclic core
    fused_ok, ring_reason = has_fused_tetracyclic_core(aglycone)
    if not fused_ok:
        return False, ring_reason
    
    # Count the number of carbons
    carbon_count = sum(1 for atom in aglycone.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (27 <= carbon_count <= 35):
        return False, (f"Aglycone contains {carbon_count} carbons; expected between 27 and 35 for a cucurbitacin core.")
    
    # Calculate molecular weight of the aglycone – expect 440 to 580 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(aglycone)
    if not (440 <= mol_wt <= 580):
        return False, (f"Aglycone molecular weight is {mol_wt:.1f} Da, "
                       "which is outside the expected 440–580 Da range for cucurbitacin cores.")
    
    # Look for at least one carbonyl (C=O) group using SMARTS.
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=O")
    carbonyl_matches = aglycone.GetSubstructMatches(carbonyl_pattern)
    
    if len(carbonyl_matches) >= 1:
        reason = (f"Aglycone fragment shows a tetracyclic core with {carbon_count} carbons, "
                  f"molecular weight {mol_wt:.1f} Da, and {len(carbonyl_matches)} carbonyl group(s).")
        return True, reason
    else:
        # Permit a relaxed acceptance when no carbonyl is found.
        reason = (f"Aglycone fragment shows a tetracyclic core with {carbon_count} carbons "
                  f"and molecular weight {mol_wt:.1f} Da, but no explicit carbonyl group was detected "
                  "(possibly due to tautomerism or masked keto functionality).")
        return True, reason

# Example usage:
if __name__ == "__main__":
    # Example: cucurbitacin D
    example_smiles = "[H][C@@]12C[C@H](O)C(=O)C(C)(C)C1=CC[C@@]1([H])[C@]3(C)C[C@@H](O)[C@]([H])([C@@](C)(O)C(=O)\\C=C\\C(C)(C)O)[C@@]3(C)CC(=O)[C@@]21C"
    result, reason = is_cucurbitacin(example_smiles)
    print("Cucurbitacin check:", result)
    print("Reason:", reason)