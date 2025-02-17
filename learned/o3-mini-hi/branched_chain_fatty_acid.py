"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: branched-chain fatty acid
A branched-chain fatty acid (BCFA) is defined as any fatty acid with one or more alkyl substituents 
on its parent hydrocarbon chain (usually a carboxylic acid with a long alkyl chain that contains at least one branch).
The following implementation:
  1. Uses a SMARTS pattern for the carboxylic acid that accepts both protonated and deprotonated forms.
  2. Finds the “R group” attached to the acid carbon.
  3. Uses a depth‐first search (DFS) to “grow” the longest carbon chain from that R group.
  4. If the chain is sufficiently long and represents a majority of the molecule’s carbon atoms 
     and if any carbon in the chain carries an extra carbon (branch) not in the longest chain, then it is classified as a BCFA.
  
Note: This approach is heuristic and may not perfectly classify every example in a difficult dataset.
"""

from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid (BCFA) based on its SMILES string.
    
    Approach:
      1. Check for a carboxylic acid moiety by means of a SMARTS pattern that matches both –COOH and –COO–.
      2. From the acid carbon, locate the R group (alkyl chain) that is attached (ignore oxygen neighbors).
      3. Find the longest connected carbon-only chain starting from that neighbor using DFS.
      4. Require that the chain is long enough and dominates the molecule’s carbon count.
      5. For atoms along that chain, check for at least one extra carbon substituent (a branch).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a branched-chain fatty acid, False otherwise.
        str: A reason explaining the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Use a SMARTS pattern for the acid that accepts both protonated and deprotonated forms.
    # [CX3](=O)[O;H1,O-] will match either -C(=O)O or -C(=O)[O-].
    acid_smarts = Chem.MolFromSmarts("[CX3](=O)[O;H1,O-]")
    acid_matches = mol.GetSubstructMatches(acid_smarts)
    if not acid_matches:
        return False, "No carboxylic acid group found; not a fatty acid."

    # For simplicity, assume the first acid match is the fatty acid head.
    # In a carboxylic acid group the first atom is the acid carbon.
    acid_idx = acid_matches[0][0]
    acid_atom = mol.GetAtomWithIdx(acid_idx)

    # Identify the R group attached to the acid carbon. 
    # It should be a carbon (atomic num 6) (ignore oxygens which are part of the acid group).
    rgroup = None
    for nbr in acid_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6:
            rgroup = nbr
            break
    if rgroup is None:
        return False, "No alkyl chain (R group) attached to the carboxyl carbon; not a fatty acid."
    
    # Create a list of all carbon atom indices in the whole molecule.
    all_carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(all_carbons) < 3:
        return False, "Too few carbon atoms to be a fatty acid."
        
    # We'll now search for the longest acyclic connected carbon chain starting from the R group.
    # (We do not exclude ring atoms so that linear unsaturated chains are not missed.)
    carbon_graph = {}
    for idx in all_carbons:
        atom = mol.GetAtomWithIdx(idx)
        # Only consider carbon neighbors.
        neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        carbon_graph[idx] = neighbors

    # DFS helper: from a starting carbon, find the longest chain.
    def dfs(current, visited):
        best_path = [current]
        for nbr in carbon_graph[current]:
            # Do not follow back into the acid carbon (if it is not our starting node) nor revisit nodes.
            if nbr in visited:
                continue
            path = [current] + dfs(nbr, visited | {nbr})
            if len(path) > len(best_path):
                best_path = path
        return best_path

    start_idx = rgroup.GetIdx()
    main_chain = dfs(start_idx, {start_idx})
    
    # For a fatty acid, the chain should be of minimal length.
    if len(main_chain) < 3:
        return False, "The fatty acyl chain appears too short."
    
    # To avoid misclassifying complex molecules (such as peptides) that happen to have a carboxyl group plus a small alkyl part,
    # require that the chain (plus acid carbon) represents a significant fraction of the molecule's carbons.
    # We consider the acid carbon as part of the acyl unit.
    acyl_chain_atoms = set(main_chain) | {acid_idx}
    ratio = len(acyl_chain_atoms) / len(all_carbons)
    if ratio < 0.5:
        return False, "The fatty acyl chain is not the dominant carbon skeleton in the molecule."

    # Now check for branch(s): go through each atom in the main chain and see if it has additional carbon neighbors
    # that are not part of the main chain. (Ignore the acid carbon if found.)
    branch_found = False
    for atom_idx in main_chain:
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            # Ignore if neighbor is in the main chain or is the acid carbon.
            if nbr_idx in main_chain or nbr_idx == acid_idx:
                continue
            branch_found = True
            break
        if branch_found:
            break

    if not branch_found:
        return False, "No alkyl substituents (branch) found on the fatty acyl chain."
    
    return True, "Contains a carboxylic acid group with a fatty acyl chain (dominant carbon skeleton) that has an alkyl branch."

# Example usage:
if __name__ == "__main__":
    # Try a few examples. For instance, isobutyric acid ("CC(C)C(O)=O") is a known BCFA.
    examples = [
        ("CC(C)C(O)=O", "isobutyric acid"),
        ("CCCC(C)C(O)=O", "2-methylvaleric acid"),
        ("OC(=O)CCC(CCCC)CC", "false positive candidate (4-Ethyloctanoic acid)"),
        ("O=C1OC(C(=O)O[C@@]23...)", "complex molecule example")
    ]
    
    for smi, name in examples:
        result, reason = is_branched_chain_fatty_acid(smi)
        print(f"SMILES: {smi}\nNAME: {name}\nResult: {result}\nReason: {reason}\n{'-'*60}")