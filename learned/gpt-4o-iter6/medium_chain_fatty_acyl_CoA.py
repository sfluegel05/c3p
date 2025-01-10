"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for thioester pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage to CoA found"

    # Broadened Coenzyme A moiety pattern focusing on adenine, ribose, phosphates, and pantetheine
    coa_moiety_pattern = Chem.MolFromSmarts(
        "NC(C=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(O)=O"
    )
    if not mol.HasSubstructMatch(coa_moiety_pattern):
        return False, "CoA moiety structure not matched"

    # Calculate the fatty acid chain length bound to the thioester group
    fatty_chain_length = calculate_chain_length(mol, thioester_pattern)

    # Medium-chain fatty acid range is about 6-12 carbons
    if not (6 <= fatty_chain_length <= 12):
        return False, f"Aliphatic chain length of {fatty_chain_length} not within medium-chain range (6-12 carbons)"
    
    return True, "Molecule is a medium-chain fatty acyl-CoA with proper CoA moiety and chain length"


def calculate_chain_length(mol, pattern):
    """
    Uses breadth-first search to calculate the longest carbon chain length starting from the thioester carbon.

    Args:
        mol: molecule object obtained from RDKit.
        pattern: SMARTS pattern that identifies thioester linkage.

    Returns:
        Longest chain length found starting at thioester carbon.
    """
    # Find the carbon linked to sulfur in thioester
    match_idx = mol.GetSubstructMatch(pattern)
    if not match_idx:
        return 0  # Failed to find starting point

    start_atom = match_idx[1]  # S carbon in C(=O)SCC

    # Use BFS to determine the longest carbon chain from start_atom
    def bfs_max_chain_length(start_idx):
        queue = [(start_idx, 0)]
        visited = set()
        max_length = 0

        while queue:
            current_idx, length = queue.pop(0)
            visited.add(current_idx)
            max_length = max(max_length, length)

            for neighbor in mol.GetAtomWithIdx(current_idx).GetNeighbors():
                n_idx = neighbor.GetIdx()
                if n_idx not in visited and neighbor.GetAtomicNum() == 6:  # Check it's a carbon
                    queue.append((n_idx, length + 1))
        
        return max_length
    
    return bfs_max_chain_length(start_atom)

# Example test
example_smiles = "CCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
print(is_medium_chain_fatty_acyl_CoA(example_smiles))