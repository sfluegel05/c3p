"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is a straight-chain fatty acid with exactly 18 carbons, two C=C double bonds,
    and a carboxylic group at one end.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Perform DFS or BFS to check if there's a straight chain of 18 carbons
    def find_longest_carbon_chain(atom):
        visited = [False] * mol.GetNumAtoms()
        max_depth = [0]

        def dfs(current, depth):
            visited[current.GetIdx()] = True
            if current.GetAtomicNum() == 6:  # Carbon atom
                max_depth[0] = max(max_depth[0], depth)
                for neighbor in current.GetNeighbors():
                    if not visited[neighbor.GetIdx()]:
                        dfs(neighbor, depth + 1)

        dfs(atom, 0)
        return max_depth[0]

    # Start BFS or DFS from each carbon atom
    max_c_chain = max(find_longest_carbon_chain(atom) for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if max_c_chain < 17:  # Because we're counting bonds, 18 carbons means 17 bonds
        return False, f"Has chain length {max_c_chain+1}, expected 18-carbon chain"
    
    # Count double bonds (consider them only within the main chain detected if possible)
    double_bond_pattern = Chem.MolFromSmarts("[C]=[C]")
    double_bonds = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bonds) != 2:
        return False, f"Found {len(double_bonds)} C=C double bonds, expected 2"

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Missing carboxylic acid group"

    # If all checks pass, it's an octadecadienoic acid
    return True, "Contains an 18-carbon chain with 2 C=C double bonds and a carboxylic acid group"

# Sample usage
smiles = "CCCCCC\C=C\C=C/CCCCCCCC(O)=O"
print(is_octadecadienoic_acid(smiles))