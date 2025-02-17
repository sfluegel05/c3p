"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: Any branched‐chain fatty acid containing methyl branches only.
The function is_methyl_branched_fatty_acid takes a SMILES string as input and returns a boolean plus a reason.
"""
from rdkit import Chem

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl–branched fatty acid (i.e., a branched–chain fatty acid
    with methyl branches only) based on its SMILES string.

    Algorithm:
      1. Parse molecule and require acyclic structure.
      2. Identify a carboxylic acid group (SMARTS: "C(=O)[O;H,-]") and use its carbon as start.
      3. Add explicit hydrogens.
      4. Build a connectivity graph of carbon atoms only.
      5. Starting at the acid carbon, compute the longest carbon-only path (candidate main chain).
      6. The chain must be at least 5 carbons long.
      7. For each carbon in the chain, examine its heavy-atom neighbors not in the main chain:
         • If the neighbor is carbon, it must form a methyl group (terminal CH3 only; i.e. having no other heavy 
           atoms besides the parent).
         • If the substituent is not carbon, then on the acid carbon allow O (as in -COOH) and on a terminal
           (non–acid) carbon allow a hydroxyl group (-OH). No other heteroatoms are permitted.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule qualifies; otherwise False.
        str: Explanation of the classification.
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that molecule is acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; not an acyclic fatty acid"
    
    # Look for a carboxylic acid group.
    acid_smarts = "C(=O)[O;H,-]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    if not acid_matches:
        return False, "No carboxylic acid group found; not a fatty acid"
    # Use first match; acid carbon is first atom in the match.
    acid_carbon = acid_matches[0][0]
    
    # Add explicit hydrogens so that we can count them accurately.
    mol_h = Chem.AddHs(mol)
    
    # Build a connectivity map over carbon atoms only (skip hydrogens)
    carbon_indices = {atom.GetIdx() for atom in mol_h.GetAtoms() if atom.GetAtomicNum() == 6}
    if acid_carbon not in carbon_indices:
        # sometimes the index may shift in hydrogens addition but we can try to map by atom symbol:
        return False, "Acid carbon not found among heavy atoms"
    
    # Build adjacency list: for each carbon, connect to neighboring carbons.
    carbon_adj = {idx: set() for idx in carbon_indices}
    for idx in carbon_indices:
        atom = mol_h.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                carbon_adj[idx].add(nbr.GetIdx())
                
    # Using DFS to find the longest simple carbon chain starting at the acid carbon.
    def longest_path(start, parent, visited):
        best = [start]
        for nbr in carbon_adj[start]:
            if nbr == parent or nbr in visited:
                continue
            path = [start] + longest_path(nbr, start, visited | {nbr})
            if len(path) > len(best):
                best = path
        return best
    
    main_chain = longest_path(acid_carbon, None, {acid_carbon})
    if len(main_chain) < 5:
        return False, "Main carbon chain is too short to be a fatty acid"
    
    # Prepare a set for fast membership check.
    main_chain_set = set(main_chain)
    
    # Helper: Check if a branch carbon (heavy atom already known to be carbon) is a methyl group.
    # A valid methyl branch should only be bonded to the backbone carbon among heavy atoms
    # and must have exactly 3 attached hydrogens.
    def is_methyl_branch(branch_idx):
        branch_atom = mol_h.GetAtomWithIdx(branch_idx)
        heavy_neighbors = [nbr for nbr in branch_atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 1:  # only the parent backbone carbon allowed
            return False
        if branch_atom.GetTotalNumHs() != 3:
            return False
        return True

    # Helper: Determine whether a hetero substituent is allowed.
    # We ignore hydrogen atoms here.
    def allowed_hetero(nbr, current_idx, terminal):
        sym = nbr.GetSymbol()
        # For acid carbon, oxygen is allowed (part of the carboxyl group).
        if current_idx == acid_carbon:
            if sym == "O":
                return True
            return False
        # For terminal backbone carbon (non–acid); allow an -OH group.
        if terminal:
            if sym == "O":
                # Check that this oxygen is -OH: i.e., bound to only one heavy atom (the backbone)
                heavy_neigh = [n for n in nbr.GetNeighbors() if n.GetAtomicNum() != 1]
                # Also, -OH should have one hydrogen.
                if len(heavy_neigh) == 1 and nbr.GetTotalNumHs() == 1:
                    return True
            return False
        # For interior carbons, no heteroatom is allowed.
        return False

    branch_count = 0
    # Now, inspect each atom in the main chain.
    for i, idx in enumerate(main_chain):
        atom = mol_h.GetAtomWithIdx(idx)
        # Determine if this is a terminal backbone carbon (but acid carbon is special)
        terminal = (i == len(main_chain) - 1) and (idx != acid_carbon)
        for nbr in atom.GetNeighbors():
            # Skip hydrogen atoms completely.
            if nbr.GetAtomicNum() == 1:
                continue
            nbr_idx = nbr.GetIdx()
            # If neighbor is part of the backbone chain, skip.
            if nbr_idx in main_chain_set:
                continue
            # If neighbor is a carbon, it must be a proper methyl branch.
            if nbr.GetAtomicNum() == 6:
                if not is_methyl_branch(nbr_idx):
                    return False, f"Branch at backbone atom {idx} is not a methyl group"
                branch_count += 1
            else:
                # For non-carbon heavy atoms, check if allowed.
                if not allowed_hetero(nbr, idx, terminal):
                    return False, f"Backbone carbon at index {idx} has disallowed substituent {nbr.GetSymbol()}"
    msg = (f"CORRECT: Methyl-branched fatty acid with main chain length {len(main_chain)} "
           f"and {branch_count} methyl branch(es)")
    return True, msg

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "CC(C)CCCCCCCCCCCCCCCCC(O)=O",           # 18-methylnonadecanoic acid
        "CC(CCCCCCC/C=C/C(=O)O)C",                # (E)-11-methyldodec-2-enoic acid
        "CC(C)CCCCCCCCCCCCCCCCCCCCCCC(O)=O",      # 24-methylpentacosanoic acid
        "OC(=O)CC(C(C)(C)C)=C",                   # 3-tert-Butyl-3-butenoic acid
        "OC(=O)CC(C)C=C",                        # 3-methyl-4-pentenoic acid
        "CCC\\C(=C/CC)C(O)=O",                    # 2-n-Propyl-2-pentenoic acid
        "CC(C)CCCCCCCCC(O)=O",                    # 10-methylundecanoic acid
        "CC(C)(C)C(O)=O",                        # pivalic acid (should be too short)
        "CC[C@@H](C)C(O)=O",                      # (R)-2-methylbutyric acid (should be too short)
        "OC(=O)CCC(C)=C",                        # 4-methyl-4-pentenoic acid
        "C(CCCCCCC(CC)C)CCCCC(O)=O",              # 13-methylpentadecanoic acid
        "CC(C)CCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # 28-methylnonacosanoic acid
        "CC(C)CCCCCCCCCCCCCCCCCCC(O)=O",          # 20-methylhenicosanoic acid
        "OC(CO)CCCCCCCCCCCCCC(O)=O",              # omega-hydroxy-15-methylpalmitic acid
        "OC(=O)\\C=C\\C(C)C",                      # 4-Methyl-2-pentenoic acid
        "OC(=O)C(CCCC(C)C)C",                      # 2,6-dimethylheptanoic acid
        "CC(C)CCCCCCCCCCCCCC(O)=O",               # isoheptadecanoic acid
        "CC(C)CCCCCC(O)=O",                       # 7-methyloctanoic acid
    ]
    
    for smi in test_smiles:
        res, reason = is_methyl_branched_fatty_acid(smi)
        print(f"{smi}  ->  {res}: {reason}")