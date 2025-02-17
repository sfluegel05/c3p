"""
Classifies: CHEBI:26333 prostaglandin
"""
#!/usr/bin/env python
"""
Classifies: Prostanoid derivatives (prostaglandins)
Definition: Naturally occurring compounds derived from the parent C20 acid, prostanoic acid.
We introduce two routes:
  A. If the molecule contains an isolated (nonfusion) cyclopentane ring (all carbons) then at least one substituent “branch”
     off the ring must be long (≥5 carbons). Optionally, the presence of an acid/ester (or amide) carbonyl (C(=O)O, C(=O)N)
     reinforces the decision.
  B. If no such ring is found, then if the molecule is mainly acyclic (or has no rings) and contains a carboxylic acid motif,
     and is long (≥20 carbons), we also classify it as prostaglandin.
If none of these conditions is met, the function returns False with an explanation.
"""

from rdkit import Chem

def is_prostaglandin(smiles: str):
    """
    Determines whether the given SMILES string is likely a prostaglandin derivative.
    The procedure applies two routes:
      Route A. Look for an isolated cyclopentane ring (all carbon atoms and non‐fused) that is decorated by at least
               one long substituent branch (≥5 carbons). Extra functional groups (acid/ester/amide carbonyls) are noted.
      Route B. In the absence of a cyclopentane core, if the molecule is largely acyclic, contains a carboxylate/acid
               motif, and has enough carbons (≥20), then classify as prostaglandin.
    Args:
      smiles (str): SMILES string of the molecule.
    Returns:
      (bool, str): Classifcation (True for prostaglandin derivative; False otherwise) and a reason.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # If molecule has multiple fragments, use the largest.
    try:
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if len(frags) > 1:
            mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    except Exception:
        pass

    # Count carbon atoms.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(carbon_atoms)
    if not (17 <= c_count <= 32):
        return False, f"Total carbon count {c_count} outside expected range for prostaglandins (17-32)"
    
    # Check for acid/ester/amide carbonyl motif.
    # We use two SMARTS patterns:
    #    acid/ester: C(=O)O
    #    amide: C(=O)N
    acid_ester = Chem.MolFromSmarts("C(=O)[O,N]")
    has_acid_like = mol.HasSubstructMatch(acid_ester)

    # Check ring information.
    ring_info = mol.GetRingInfo().AtomRings()
    candidate_rings = []  # to hold isolated cyclopentane rings (non-fused, all C)
    # Count ring membership per atom.
    ring_membership = {atom.GetIdx(): 0 for atom in mol.GetAtoms()}
    for ring in ring_info:
        for idx in ring:
            ring_membership[idx] += 1
            
    for ring in ring_info:
        if len(ring) == 5:
            # Check that all atoms are carbon.
            if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                # Check that every atom in the ring belongs to exactly one ring (non-fused)
                if all(ring_membership[idx] == 1 for idx in ring):
                    candidate_rings.append(set(ring))
                    
    # Helper function to recursively collect branch atoms starting from a given atom (by idx)
    # that are not part of the specified core.
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

    # Route A: if a candidate cyclopentane ring exists, check for a long branch.
    for core in candidate_rings:
        # For atoms in the ring, collect neighbors not in the ring.
        substituent_idxs = set()
        for idx in core:
            atom = mol.GetAtomWithIdx(idx)
            for nb in atom.GetNeighbors():
                if nb.GetIdx() not in core:
                    substituent_idxs.add(nb.GetIdx())
        # For each branch, measure its carbon count.
        long_branch_found = False
        for sub_idx in substituent_idxs:
            try:
                branch = get_branch(sub_idx, core)
            except Exception as e:
                return False, f"Error processing a substituent branch: {e}"
            branch_carbon_count = sum(1 for idx in branch if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if branch_carbon_count >= 5:
                long_branch_found = True
                break
        if long_branch_found:
            reason = "Contains an isolated cyclopentane core with a long substituent branch"
            if has_acid_like:
                reason += " and an acid/ester/amide carbonyl motif."
            else:
                reason += " (no acid/ester/amide motif detected, so likely an alcohol derivative)."
            return True, reason
    
    # Route B: In case no eligible cyclopentane ring is found.
    # Check if the molecule is largely acyclic (or has no small rings) and exhibits an acid functionality.
    if not ring_info or all(len(ring) > 6 for ring in ring_info):
        acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-]")  # carboxyl/function motif
        if mol.HasSubstructMatch(acid_pattern) and c_count >= 20:
            return True, ("Acyclic (or non-cyclized) long-chain molecule with a carboxylate/acid motif; "
                          "possible prostaglandin derivative (e.g. nonacosanoic acid).")
    
    return False, "No prostaglandin-like core was found: missing isolated cyclopentane with long branch or acyclic acid motif."


# Example usage:
if __name__ == "__main__":
    test_examples = {
        "prostaglandin A2": "CCCC[C@H](O)\\C=C\\[C@H]1C=CC(=O)[C@@H]1C\\C=C/CCCC(O)=O",
        "nonacosanoic acid": "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # Acyclic acid
        "prostaglandin E3": "CC\\C=C/C[C@H](O)\\C=C\\[C@H]1[C@H](O)CC(=O)[C@@H]1C\\C=C/CCCC(O)=O",
        "prostaglandin F1alpha alcohol": "C(O)CCCCCC[C@@H]1[C@H]([C@@H](C[C@@H]1O)O)/C=C/[C@H](CCCCC)O",
        "prostaglandin F2alpha dimethylamine": "CCCC[C@H](O)\\C=C\\[C@H]1[C@H](O)C[C@H](O)[C@@H]1C\\C=C/CCCCN(C)C",
    }
    for name, s in test_examples.items():
        result, reason = is_prostaglandin(s)
        print(f"Name: {name}\nSMILES: {s}\nResult: {result}\nReason: {reason}\n")