"""
Classifies: CHEBI:26333 prostaglandin
"""
#!/usr/bin/env python
"""
Classifies: Naturally occurring prostaglandin derivatives
Definition: Compounds derived from the parent C20 acid, prostanoic acid.
We now require that either (A) the molecule contains a 5‐membered (prostan) core 
– where the core is either isolated or fused to at most one additional ring – with at least one 
long substituent chain (≥5 carbons) and optionally showing an acid/ester/amide carbonyl motif, or 
(B) if no such cyclic core is present then (if the molecule is mainly acyclic, contains a carboxylic acid group, 
and has a total carbon count in the expected range of the C20 parent, i.e. 17–23 carbons) it is classified as a prostaglandin derivative.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines whether the given SMILES string is likely a prostaglandin derivative.
    Two routes are applied:
    
    Route A: Look for at least one five-membered ring (all carbon) 
             that is not highly fused (each atom participates in at most 2 rings)
             and that is decorated by at least one substituent branch (chain) that has ≥5 carbon atoms.
             Optionally, if a carbonyl motif (acid, ester or amide: C(=O)[O,N]) is present in the molecule,
             that reinforces the decision.
    
    Route B: If no such cyclic prostaglandin core is found, then if the molecule is largely acyclic
             (or only contains rings of size ≥7) and its total number of carbons is in the 17–23 range 
             (i.e. near the parent C20 skeleton) and it contains a carboxylic acid motif,
             then it is classified as prostaglandin.
             
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): A tuple with a boolean (True if the molecule is classified as prostaglandin)
                     and a reason string.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # If molecule has multiple fragments, work on the largest fragment.
    try:
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if len(frags) > 1:
            mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    except Exception:
        pass

    # Count total carbon atoms.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(carbons)
    # For route B, we require the carbon count to be near that of prostanoic acid (C20), say 17-23.
    in_expected_range = (17 <= c_count <= 23)

    # SMARTS for acid/ester/amide carbonyl group, used as reinforcing evidence.
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)[O,N]")
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)
    # SMARTS for a carboxylic acid moiety.
    acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")

    # Get ring information.
    ring_info = mol.GetRingInfo().AtomRings()
    atom_ring_membership = {atom.GetIdx(): 0 for atom in mol.GetAtoms()}
    for ring in ring_info:
        for idx in ring:
            atom_ring_membership[idx] += 1

    # Find candidate 5-membered rings (all carbon) that are “not too fused” 
    # (we allow an atom to be in at most 2 rings).
    candidate_rings = []
    for ring in ring_info:
        if len(ring) == 5:
            if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                if all(atom_ring_membership[idx] <= 2 for idx in ring):
                    candidate_rings.append(set(ring))
    
    # Helper: recursively obtain the connected branch (atoms not in the core) 
    def get_branch(atom_idx, core, visited=None):
        if visited is None:
            visited = set()
        visited.add(atom_idx)
        branch_atoms = {atom_idx}
        atom = mol.GetAtomWithIdx(atom_idx)
        for nb in atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx in core or nb_idx in visited:
                continue
            branch_atoms.update(get_branch(nb_idx, core, visited))
        return branch_atoms

    # Route A: Look for a candidate cyclopentane (or five‐membered) core decorated with a long branch.
    for core in candidate_rings:
        # Find substituents: atoms adjacent to the ring that are not in the ring.
        substituents = set()
        for idx in core:
            atom = mol.GetAtomWithIdx(idx)
            for nb in atom.GetNeighbors():
                if nb.GetIdx() not in core:
                    substituents.add(nb.GetIdx())
        # Check each substituent branch for long carbon chain (≥5 carbon atoms).
        for sub_idx in substituents:
            try:
                branch_indices = get_branch(sub_idx, core)
            except Exception as e:
                return False, f"Error processing substituent branch: {e}"
            branch_carbons = sum(1 for i in branch_indices if mol.GetAtomWithIdx(i).GetAtomicNum() == 6)
            if branch_carbons >= 5:
                reason = "Contains a 5‐membered prostaglandin core with a long substituent branch"
                if has_carbonyl:
                    reason += " and an acid/ester/amide carbonyl motif."
                else:
                    reason += " (no acid/ester/amide motif detected)."
                return True, reason

    # Route B: For acyclic or non-prostan molecules.
    # We allow only molecules that are either completely acyclic or whose rings are all large (≥7 members).
    if (not ring_info or all(len(ring) >= 7 for ring in ring_info)) and in_expected_range and mol.HasSubstructMatch(acid_pattern):
        reason = ("Acyclic (or non-cyclized) molecule with total carbon count of {} and a carboxylic acid motif; "
                  "possible prostaglandin derivative.").format(c_count)
        return True, reason

    return False, "No prostaglandin-like core was found: missing a suitable five-membered ring with a long branch or the acyclic acid criteria."

# Example usage:
if __name__ == "__main__":
    test_examples = {
        "prostaglandin A2": "CCCC[C@H](O)\\C=C\\[C@H]1C=CC(=O)[C@@H]1C\\C=C/CCCC(O)=O",
        "nonacosanoic acid": "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",
        "prostaglandin E3": "CC\\C=C/C[C@H](O)\\C=C\\[C@H]1[C@H](O)CC(=O)[C@@H]1C\\C=C/CCCC(O)=O",
        "6alpha-Prostaglandin I1": "[C@@H]1([C@@H]([C@@H]2[C@H](C1)O[C@@H](C2)CCCCC(O)=O)/C=C/[C@H](CCCCC)O)O",
        "prostaglandin F1alpha alcohol": "C(O)CCCCCC[C@@H]1[C@H]([C@@H](C[C@@H]1O)O)/C=C/[C@H](CCCCC)O",
        "prostaglandin H3": "[C@@H]1([C@@H]2C[C@H]([C@@H]1/C=C/[C@H](C/C=C\\CC)O)OO2)C/C=C\\CCCC(O)=O",
    }
    for name, s in test_examples.items():
        result, reason = is_prostaglandin(s)
        print(f"Name: {name}\nSMILES: {s}\nResult: {result}\nReason: {reason}\n")