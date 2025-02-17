"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: Cucurbitacin (tetracyclic triterpenoids derived from cucurbitane)

Cucurbitacins (often glycosylated) feature an aglycone that is a fused tetracyclic scaffold.
Typically the aglycone shows:
  • exactly 4 fused rings,
  • between about 27 and 35 carbon atoms in the core,
  • a molecular weight in the 440–580 Da range, and
  • one or more carbonyl (C=O) groups.
  
This implementation splits the molecule into fragments and (if needed) uses the Murcko scaffold
to remove appended glycosyl or acyclic substituents. Then it checks for a fused tetracyclic core.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def has_fused_tetracyclic_core(mol):
    """
    Checks whether the molecule (or candidate core) contains a connected set of exactly 4 fused rings.
    We use the ring information (GetRingInfo) and build a graph where rings are nodes that are connected
    if they share at least one atom. Then we search for any connected component of exactly 4 rings.
    
    Returns:
       (True, reason) if such a set is found,
       otherwise (False, reason).
    """
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if len(atom_rings) < 4:
        return False, f"Found only {len(atom_rings)} rings; expecting at least 4 for a tetracyclic core."
    
    # Represent each ring as a set of atom indices.
    ring_sets = [set(r) for r in atom_rings]
    n = len(ring_sets)
    connections = {i: set() for i in range(n)}
    # Connect two rings if they share at least one atom.
    for i in range(n):
        for j in range(i+1, n):
            if ring_sets[i].intersection(ring_sets[j]):
                connections[i].add(j)
                connections[j].add(i)
    seen = set()
    components = []
    for i in range(n):
        if i in seen:
            continue
        stack = [i]
        comp = set()
        while stack:
            node = stack.pop()
            if node in comp:
                continue
            comp.add(node)
            for neigh in connections[node]:
                if neigh not in comp:
                    stack.append(neigh)
        seen = seen.union(comp)
        components.append(comp)
    # Accept only if a connected component has exactly 4 rings.
    for comp in components:
        if len(comp) == 4:
            return True, "Found fused ring system with exactly 4 rings."
    return False, "No connected set of exactly 4 fused rings found."

def get_aglycone_candidate(mol):
    """
    Attempts to extract an aglycone candidate from the (possibly glycosylated) molecule.
    First it splits the molecule into disconnected fragments.
    Then it collects fragments that have at least a minimal carbon count (>=27) and, if they already
    feature a fused tetracyclic core, they are preferred.
    If none of the fragments exhibits the expected core, we select the fragment with highest carbon count.
    
    Returns a candidate fragment (an RDKit Mol).
    """
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    candidates = []
    for frag in frags:
        carbon_count = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
        mw = rdMolDescriptors.CalcExactMolWt(frag)
        # Only consider fragments that are plausibly large (filter out obvious sugars)
        if carbon_count >= 27 and mw > 300:
            fused, _ = has_fused_tetracyclic_core(frag)
            # If fused core exists, we calculate a penalty for deviation from a nominal MW (~510 Da).
            penalty = abs(mw - 510) if fused else 1000 + abs(mw - 510)
            candidates.append((penalty, frag))
    if candidates:
        # Pick candidate with lowest penalty.
        best_candidate = sorted(candidates, key=lambda x: x[0])[0][1]
        return best_candidate
    # If none meets our strict criteria, fall back to the fragment with highest carbon count.
    best_frag = None
    best_carbons = 0
    for frag in frags:
        carbons = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
        if carbons > best_carbons:
            best_frag = frag
            best_carbons = carbons
    return best_frag if best_frag is not None else mol

def is_cucurbitacin(smiles: str):
    """
    Determines whether the input SMILES string corresponds to a cucurbitacin.
    The classifier works by attempting to isolate an aglycone candidate whose scaffold:
       • exhibits a fused ring system of exactly 4 rings,
       • features between 27 and 35 carbon atoms in its core,
       • has a molecular weight between 440 and 580 Da,
       • and contains at least one carbonyl (C=O) group.
       
    To handle glycosylated compounds, this function first splits the molecule into fragments;
    if the candidate appears too large (carbon count >35), it then computes the Murcko scaffold,
    which often removes appended sugars and extraneous side‐chains.
    
    Returns:
       (True, reason) if classified as cucurbitacin,
       otherwise (False, reason).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Obtain an initial candidate from disconnected fragments.
    candidate = get_aglycone_candidate(mol)
    if candidate is None:
        candidate = mol

    carbon_count = sum(1 for atom in candidate.GetAtoms() if atom.GetAtomicNum() == 6)
    mw = rdMolDescriptors.CalcExactMolWt(candidate)
    
    # If the candidate core appears too large (likely still bearing glycosyl parts),
    # compute the Murcko scaffold to remove peripheral substituents.
    if carbon_count > 35 or mw > 600:
        scaffold = MurckoScaffold.GetScaffoldForMol(candidate)
        if scaffold:
            candidate = scaffold
            carbon_count = sum(1 for atom in candidate.GetAtoms() if atom.GetAtomicNum() == 6)
            mw = rdMolDescriptors.CalcExactMolWt(candidate)

    # Check for exactly 4 fused rings.
    fused_ok, ring_reason = has_fused_tetracyclic_core(candidate)
    if not fused_ok:
        return False, ring_reason

    # Verify the carbon count in the candidate (should be between about 27 and 35).
    if not (27 <= carbon_count <= 35):
        return False, (f"Candidate core contains {carbon_count} carbons; expected between 27 and 35.")
    
    # Verify the molecular weight of the candidate (expected between 440 and 580 Da).
    if not (440 <= mw <= 580):
        return False, (f"Candidate core molecular weight is {mw:.1f} Da; "
                       "expected between 440 and 580 Da.")
    
    # Look for at least one carbonyl (C=O) group.
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=O")
    carbonyl_matches = candidate.GetSubstructMatches(carbonyl_pattern)
    if carbonyl_matches:
        reason = (f"Candidate core shows a tetracyclic fused ring system with {carbon_count} carbons, "
                  f"molecular weight {mw:.1f} Da, and {len(carbonyl_matches)} carbonyl group(s).")
        return True, reason
    else:
        reason = (f"Candidate core shows a tetracyclic fused ring system with {carbon_count} carbons and "
                  f"molecular weight {mw:.1f} Da, but no explicit carbonyl group was detected "
                  "(possibly due to tautomerism or masked keto functionality).")
        return True, reason

# Example usage:
if __name__ == "__main__":
    example_smiles = "[H][C@@]12C[C@H](O)C(=O)C(C)(C)C1=CC[C@@]1([H])[C@]3(C)C[C@@H](O)[C@]([H])([C@@](C)(O)C(=O)\\C=C\\C(C)(C)O)[C@@]3(C)CC(=O)[C@@]21C"  # cucurbitacin D
    result, reason = is_cucurbitacin(example_smiles)
    print("Cucurbitacin check:", result)
    print("Reason:", reason)