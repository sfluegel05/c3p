"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Classifies a molecule as a very long-chain fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA, False otherwise.
        str: Reason for classification.
    """
    # Convert the SMILES string to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Substructure pattern for CoA (coenzyme A) part
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA structure found"

    # Substructure pattern for the carbonyl group linking the CoA and fatty acyl chain
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)SC")
    carbonyl_match = mol.GetSubstructMatch(carbonyl_pattern)

    if not carbonyl_match:
        return False, "No carbonyl linked to CoA structure found"

    # Starting from carbonyl carbon, find longest consecutive chain of carbon atoms
    carbonyl_carbon = carbonyl_match[0]  # carbon of the carbonyl group

    visited = set()
    longest_chain_length = 0

    def dfs(atom, length):
        nonlocal longest_chain_length
        visited.add(atom.GetIdx())
        
        # Update longest chain length
        if length > longest_chain_length:
            longest_chain_length = length

        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                dfs(neighbor, length + 1)

    # Begin DFS from carbon atom adjacent to our identified starting point
    for atom in mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors():
        if atom.GetAtomicNum() == 6:
            dfs(atom, 1)

    if longest_chain_length <= 22:
        return False, f"Longest carbon chain is {longest_chain_length}, must be greater than 22"

    return True, "Contains CoA structure with a fatty acyl chain longer than C22"