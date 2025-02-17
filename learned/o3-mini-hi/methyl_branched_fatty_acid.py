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
    Determines if a molecule is a methyl–branched fatty acid (a branched–chain fatty acid with methyl branches only)
    based on its SMILES string.

    The algorithm:
      1. Parse the molecule and add explicit hydrogens.
      2. Check that the molecule is acyclic.
      3. Identify a carboxylic acid group (using a SMARTS that matches -C(=O)[O;H,-]).
      4. Build a carbon-only connectivity (graph) from the explicit–H molecule.
      5. Starting from the acid carbon, compute the longest possible simple path (the candidate main chain).
      6. Verify that the candidate chain is long enough (we require at least 5 carbons).
      7. For each carbon on the chain, check that any extra carbon (branch) is a terminal CH3.
         Also check that non‐carbon substituents are only allowed:
           - On the acid carbon: oxygen is allowed.
           - On the terminal (non–acid) carbon: an –OH or a quaternary nitrogen (with no H’s) is allowed.
           - On interior carbons: no heteroatom is allowed.
      8. If all conditions are met, report True with a message.

    Args:
      smiles (str): The SMILES string of the molecule.

    Returns:
      bool: True if the molecule is classified as a methyl–branched fatty acid; otherwise False.
      str: Explanation/reason for the classification.
    """
    # Parse molecule and check validity.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Require acyclic molecules.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; not an acyclic fatty acid"
    
    # Look for a carboxylic acid group.
    # This SMARTS matches protonated or deprotonated -COOH groups.
    acid_smarts = "C(=O)[O;H,-]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    if not acid_matches:
        return False, "No carboxylic acid group found; not a fatty acid"
    # Use the first match; assume the acid carbon (index 0) is the start of the main chain.
    acid_carbon = acid_matches[0][0]
    
    # Add explicit hydrogens so that we can check hydrogen counts.
    mol_h = Chem.AddHs(mol)
    
    # Build carbon-only connectivity for atoms in the explicit molecule.
    carbon_indices = {atom.GetIdx() for atom in mol_h.GetAtoms() if atom.GetAtomicNum() == 6}
    if acid_carbon not in carbon_indices:
        return False, "Acid carbon not found among carbons"
    
    # Build an adjacency list for carbon atoms.
    carbon_adj = {idx: set() for idx in carbon_indices}
    for idx in carbon_indices:
        atom = mol_h.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                carbon_adj[idx].add(nbr.GetIdx())
    
    # Helper function: recursively compute the longest path from a starting carbon.
    def get_longest_path(current, parent):
        best_path = [current]
        for nbr in carbon_adj[current]:
            if nbr == parent:
                continue
            candidate = [current] + get_longest_path(nbr, current)
            if len(candidate) > len(best_path):
                best_path = candidate
        return best_path

    main_chain = get_longest_path(acid_carbon, None)
    if len(main_chain) < 5:
        return False, "Main carbon chain is too short to be a fatty acid"
    
    # A helper to check if a branch carbon is an acceptable methyl substituent.
    def is_methyl_branch(branch_idx, parent_idx):
        # It must be a carbon that is terminal in the carbon graph.
        if len(carbon_adj[branch_idx]) != 1:
            return False
        branch_atom = mol_h.GetAtomWithIdx(branch_idx)
        # The branch must be CH3 (3 hydrogens).
        if branch_atom.GetTotalNumHs() != 3:
            return False
        return True

    # Helper to check heteroatom substituents on a backbone carbon.
    def allowed_hetero(substituent, backbone_idx, terminal=False):
        sym = substituent.GetSymbol()
        # On the acid carbon, allow oxygen (as in the -COOH group).
        if backbone_idx == acid_carbon:
            if sym == 'O':
                return True
            return False
        # On a terminal carbon (non acid), allow an -OH group.
        if terminal:
            if sym == 'O':
                # Expect it to be hydroxyl (one hydrogen only).
                num_H = sum(1 for n in substituent.GetNeighbors() if n.GetSymbol() == 'H')
                if num_H == 1:
                    return True
            # Also allow quaternary nitrogen (no hydrogens) if needed.
            if sym == 'N' and substituent.GetTotalNumHs() == 0:
                return True
            return False
        # For interior backbone carbons, do not allow heteroatoms.
        return False

    # Now, inspect each carbon on the main chain.
    chain_set = set(main_chain)
    branch_count = 0
    for i, idx in enumerate(main_chain):
        atom = mol_h.GetAtomWithIdx(idx)
        # Determine if this is a terminal carbon (but note that acid carbon is special).
        is_terminal = (i == len(main_chain)-1) and (idx != acid_carbon)
        # Iterate over all neighbors of this carbon.
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr.GetAtomicNum() == 6:
                # If neighbor is carbon not in the main chain, it is a branch.
                if nbr_idx not in chain_set:
                    if not is_methyl_branch(nbr_idx, idx):
                        return False, f"Branch at backbone atom {idx} is not a methyl group"
                    branch_count += 1
            else:
                # For non-carbon substituents, check allowed patterns.
                # Skip if the neighbor is part of the carboxylic acid already attached at acid carbon.
                # On acid carbon, any oxygen is allowed.
                if idx == acid_carbon:
                    if nbr.GetSymbol() == 'O':
                        continue
                    else:
                        return False, f"Acid carbon has disallowed substituent {nbr.GetSymbol()}"
                else:
                    if not allowed_hetero(nbr, idx, terminal=is_terminal):
                        return False, f"Backbone carbon at index {idx} has disallowed hetero substituent {nbr.GetSymbol()}"
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
        print(smi, res, reason)