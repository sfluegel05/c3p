"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: very long-chain fatty acid (chain length > C22; >C27 = ultra-long-chain)
Heuristic: must contain a carboxylic acid group and an unbroken chain of carbons 
with a total chain length (starting from the acid carbon) greater than 22.
"""
from rdkit import Chem

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid is defined as having a fatty acid chain (starting from the carboxylic acid carbon)
    with more than 22 carbons in series. (Fatty acids with more than 27 chain carbons are often called ultra-long-chain.)

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule qualifies as very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First check for presence of the carboxylic acid group.
    # We use two SMARTS patterns to capture both protonated and deprotonated forms.
    acid_smarts1 = Chem.MolFromSmarts("C(=O)[OH]")   # protonated acid
    acid_smarts2 = Chem.MolFromSmarts("C(=O)[O-]")    # deprotonated acid

    matches = mol.GetSubstructMatches(acid_smarts1) + mol.GetSubstructMatches(acid_smarts2)
    if not matches:
        return False, "No carboxylic acid group found"

    # Use the first match; the first atom in the match is the acid carbon.
    acid_carbon_idx = matches[0][0]

    # Build a graph that contains only carbon atoms.
    # We'll represent the graph as a dictionary mapping each carbon atom index (that is in the molecule)
    # to a list of adjacent carbon atom indices.
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    graph = {idx: [] for idx in carbon_idxs}
    for bond in mol.GetBonds():
        a = bond.GetBeginAtom()
        b = bond.GetEndAtom()
        if a.GetAtomicNum() == 6 and b.GetAtomicNum() == 6:
            graph[a.GetIdx()].append(b.GetIdx())
            graph[b.GetIdx()].append(a.GetIdx())

    if acid_carbon_idx not in graph:
        return False, "Acid carbon not found in the carbon graph, unexpected structure"

    # Define a recursive DFS routine to find the longest simple path (no repeated atoms)
    def dfs(node, visited):
        # Count current node in the path length.
        max_length = 1
        for neighbor in graph[node]:
            if neighbor not in visited:
                # Extend the visited set and search further
                length = 1 + dfs(neighbor, visited | {neighbor})
                if length > max_length:
                    max_length = length
        return max_length

    # Compute the longest continuous chain starting from the acid carbon.
    longest_chain_length = dfs(acid_carbon_idx, {acid_carbon_idx})

    # Check if the chain length meets the threshold (must be > 22 carbons)
    if longest_chain_length <= 22:
        return False, f"Longest carbon chain from the acid carbon has {longest_chain_length} carbons, which is not >22"
    else:
        reason = f"Longest carbon chain from the acid carbon has {longest_chain_length} carbons; qualifies as very long-chain fatty acid"
        if longest_chain_length > 27:
            reason += " (ultra-long-chain fatty acid)"
        return True, reason