"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.
    A very long-chain fatty acyl-CoA has a fatty acyl group with a chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for coenzyme A substructure
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A backbone found"

    # Find all carbon chains in the molecule
    def find_longest_chain(mol):
        # A helper function to look for the longest path of carbon atoms
        max_length = 0

        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'C':
                # Perform BFS or DFS to find the longest carbon chain starting from this atom
                stack = [(atom, 0, {atom.GetIdx()})]  # (current atom, current length, visited set)
                while stack:
                    current_atom, length, visited = stack.pop()
                    max_length = max(max_length, length)

                    for neighbor in current_atom.GetNeighbors():
                        if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() not in visited:
                            stack.append((neighbor, length + 1, visited | {neighbor.GetIdx()}))

        return max_length

    longest_chain_length = find_longest_chain(mol)

    # Check if the longest carbon chain > 22
    if longest_chain_length > 22:
        return True, f"Contains CoA backbone and fatty acyl chain length is {longest_chain_length}, which is greater than C22"
    
    return False, f"Longest carbon chain length is {longest_chain_length}, not greater than C22"