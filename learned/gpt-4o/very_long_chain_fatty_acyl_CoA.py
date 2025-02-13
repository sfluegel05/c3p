"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA 
    based on its SMILES string. Very long-chain fatty acyl-CoAs 
    have fatty acyl groups with chain lengths greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acyl-CoA,
              False otherwise
        str: Reason for classification
    """

    # Parse SMILES into RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify Coenzyme A part by locating substructures
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A (CoA) core structure found"
    
    # Locate the acyl chain part, i.e., part before thioester linkage
    # Use S atom in thioester linkage (closest to CoA pattern) to identify acyl chain
    thioester_sulfur = None
    for match in mol.GetSubstructMatches(coa_pattern):
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 16 and any(n.GetAtomicNum() == 6 for n in atom.GetNeighbors()):
                thioester_sulfur = atom
                break
    
    if thioester_sulfur is None:
        return False, "No thioester linkage found"

    # Use DFS to find longest acyl chain starting from carbon neighbors of thioester sulfur
    def find_longest_linear_chain(atom, visited=None):
        if visited is None:
            visited = set()
        max_length = 0
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum() == 6:
                visited.add(neighbor.GetIdx())
                chain_length = 1 + find_longest_linear_chain(neighbor, visited)
                max_length = max(max_length, chain_length)
                visited.remove(neighbor.GetIdx())
        return max_length

    # Start counting carbons from the carbon directly bonded to sulfur, do not include sulfur itself
    max_chain_length = 0
    for neighbor in thioester_sulfur.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:
            length = find_longest_linear_chain(neighbor)
            max_chain_length = max(max_chain_length, length)
    
    # Check the length of the identified acyl chain
    if max_chain_length > 22:
        return True, f"Fatty acyl chain length is {max_chain_length} (> 22 carbons, valid)"
    
    return False, f"Fatty acyl chain length is {max_chain_length} (<= 22 carbons, not valid)"