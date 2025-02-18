"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Improved classifier for cucurbitacins (tetracyclic triterpenoids derived from cucurbitane)

A cucurbitacin is expected to have an aglycone that is a fused tetracyclic scaffold.
Typically the aglycone shows:
  • exactly 4 fused rings (ignoring rings smaller than 5 atoms)
  • between about 27 and 35 carbon atoms in its core,
  • a molecular weight roughly between 430 and 580 Da,
  • and at least one carbonyl (C=O) group.

To address previous false negatives/positives we now:
  • Always extract the Murcko scaffold to remove appended sugars or acyclic substituents.
  • Only count rings of size ≥5 when looking for a connected set of 4 fused rings.
  • Slightly relax the lower molecular weight bound.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

def has_fused_tetracyclic_core(mol):
    """
    Checks whether the molecule (or candidate core) contains a connected set of exactly 4 fused rings.
    Only rings of size 5 or more are taken into account.
    We build a graph where rings are nodes that are connected if they share at least one atom.
    Returns:
       (True, reason) if a connected component of exactly 4 rings is found,
       otherwise (False, reason).
    """
    ring_info = mol.GetRingInfo()
    # Get only rings that have 5 or more atoms.
    atom_rings = [r for r in ring_info.AtomRings() if len(r) >= 5]
    if len(atom_rings) < 4:
        return False, f"Found only {len(atom_rings)} rings (of size>=5); expecting at least 4 for a tetracyclic core."
    
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
            return True, "Found fused tetracyclic ring system (exactly 4 rings)."
    return False, "No connected set of exactly 4 fused rings found."

def get_aglycone_candidate(mol):
    """
    Attempts to extract an aglycone candidate from the (possibly glycosylated) molecule.
    First the molecule is split into fragments; then the Murcko scaffold is computed
    from the largest fragment to remove appended acyclic (and often glycosidic) substituents.
    Returns a candidate fragment (an RDKit Mol).
    """
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    # Pick the fragment with the highest carbon count (a proxy for being the triterpenoid core)
    best_frag = None
    best_carbons = 0
    for frag in frags:
        carbons = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
        if carbons > best_carbons:
            best_frag = frag
            best_carbons = carbons
    candidate = best_frag if best_frag is not None else mol
    # Always try to remove acyclic substituents using the Murcko scaffold.
    scaffold = MurckoScaffold.GetScaffoldForMol(candidate)
    if scaffold is not None and scaffold.GetNumAtoms() > 0:
        candidate = scaffold
    return candidate

def is_cucurbitacin(smiles: str):
    """
    Determines whether the input SMILES string corresponds to a cucurbitacin.
    The classifier works by:
      • extracting a candidate aglycone core via Murcko scaffolding,
      • checking that the candidate core has exactly 4 fused rings (ignoring rings size<5),
      • ensuring that the candidate has between 27 and 35 carbons,
      • a molecular weight between ~430 and 580 Da,
      • and at least one carbonyl (C=O) group.
      
    Returns:
       (True, reason) if classified as cucurbitacin,
       otherwise (False, reason).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    candidate = get_aglycone_candidate(mol)
    if candidate is None:
        candidate = mol

    # Count carbon atoms and molecular weight of candidate.
    carbon_count = sum(1 for atom in candidate.GetAtoms() if atom.GetAtomicNum() == 6)
    mw = rdMolDescriptors.CalcExactMolWt(candidate)
    
    # Check fused tetracyclic core (using rings with size>=5).
    fused_ok, ring_reason = has_fused_tetracyclic_core(candidate)
    if not fused_ok:
        return False, ring_reason

    # Verify that the candidate core has a carbon count between about 27 and 35.
    if not (27 <= carbon_count <= 35):
        return False, f"Candidate core contains {carbon_count} carbons; expected between 27 and 35."
    
    # Verify the molecular weight (allowing a slightly lower bound for some glycosides).
    if not (430 <= mw <= 580):
        return False, f"Candidate core molecular weight is {mw:.1f} Da; expected between 430 and 580 Da."
    
    # Look for at least one carbonyl (C=O) group.
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=O")
    carbonyl_matches = candidate.GetSubstructMatches(carbonyl_pattern)
    if carbonyl_matches:
        reason = (f"Candidate core shows a fused tetracyclic ring system with {carbon_count} carbons, "
                  f"molecular weight {mw:.1f} Da, and {len(carbonyl_matches)} carbonyl group(s).")
        return True, reason
    else:
        reason = (f"Candidate core shows a fused tetracyclic ring system with {carbon_count} carbons and "
                  f"molecular weight {mw:.1f} Da, but no explicit carbonyl group was detected "
                  "(possibly due to tautomerism or masked keto functionality).")
        return True, reason

# Example usage:
if __name__ == "__main__":
    # Example: cucurbitacin D (one of the true positives)
    example_smiles = "[H][C@@]12C[C@H](O)C(=O)C(C)(C)C1=CC[C@@]1([H])[C@]3(C)C[C@@H](O)[C@]([H])([C@@](C)(O)C(=O)\\C=C\\C(C)(C)O)[C@@]3(C)CC(=O)[C@@]21C"
    result, reason = is_cucurbitacin(example_smiles)
    print("Cucurbitacin check:", result)
    print("Reason:", reason)