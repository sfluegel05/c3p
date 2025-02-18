"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:XXXXX long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A long-chain fatty acyl-CoA is a thioester formed between coenzyme A and a long-chain fatty acid (C13-C22).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find thioester groups (S-C(=O))
    thioester_pattern = Chem.MolFromSmarts('[S]-[C]=[O]')
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found"
    
    # Check for pantetheine fragment (S-C-C-N-C(=O)) indicative of CoA
    pantetheine_pattern = Chem.MolFromSmarts('S-C-C-N-C(=O)')
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Pantetheine fragment not found (not CoA?)"
    
    # Check each thioester group for chain length
    for match in thioester_matches:
        s_idx = match[0]
        c_idx = match[1]
        sulfur = mol.GetAtomWithIdx(s_idx)
        carbonyl_carbon = mol.GetAtomWithIdx(c_idx)
        
        # Find R group atom (connected to carbonyl, not S or O)
        r_group_atom = None
        for neighbor in carbonyl_carbon.GetNeighbors():
            if neighbor.GetIdx() == s_idx:
                continue
            if neighbor.GetAtomicNum() == 8:  # O in carbonyl
                continue
            r_group_atom = neighbor
            break
        if not r_group_atom:
            continue  # shouldn't happen if thioester exists
        
        # Traverse R group to count carbons
        visited = set()
        stack = [r_group_atom]
        carbon_count = 0
        while stack:
            atom = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            if atom.GetAtomicNum() != 6:  # count only carbons
                continue
            carbon_count += 1
            for neighbor in atom.GetNeighbors():
                # Avoid going back to carbonyl or into CoA
                if neighbor.GetIdx() != c_idx and neighbor.GetAtomicNum() != 16:  # 16 is S
                    stack.append(neighbor)
        
        # Check chain length (R group carbons + carbonyl carbon? Or just R group?)
        # Original fatty acid is C13-C22, so R group (excluding carbonyl) should be 12-21 carbons
        if 12 <= carbon_count <= 21:
            return True, f"Long-chain acyl group ({carbon_count + 1} carbons) attached via thioester to CoA"
    
    return False, "Acyl chain length not in C13-C22 range"