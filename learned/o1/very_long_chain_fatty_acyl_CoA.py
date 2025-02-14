"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
"""
Classifies: very long-chain fatty acyl-CoA
Definition: A fatty acyl-CoA in which the fatty acyl group has a chain length greater than C22.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA based on its SMILES string.
    A very long-chain fatty acyl-CoA is a fatty acyl-CoA molecule where the acyl chain has more than 22 carbons.
    The acyl chain should be linear (unbranched) and aliphatic (no rings).

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

    # Find the thioester bond: C(=O)-S-
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    matches = mol.GetSubstructMatches(thioester_pattern)
    if not matches:
        return False, "No thioester linkage found"
    
    # Assume the first match is the acyl-CoA linkage
    match = matches[0]
    carbonyl_c_idx = match[0]
    sulfur_idx = match[2]  # Index of the sulfur atom (should be connected to CoA)
    
    # Break the thioester bond to separate acyl chain from CoA
    cs_bond = mol.GetBondBetweenAtoms(carbonyl_c_idx, sulfur_idx)
    if cs_bond is None:
        return False, "Failed to identify thioester bond"
    
    bond_idx = cs_bond.GetIdx()
    
    # Break the thioester bond
    mol_frag = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=False)
    frags = Chem.GetMolFrags(mol_frag, asMols=True, sanitizeFrags=True)
    
    # Identify the acyl chain fragment
    acyl_chain = None
    for frag in frags:
        # Check if fragment contains the carbonyl carbon (part of acyl chain)
        if frag.HasSubstructMatch(Chem.MolFromSmarts("C(=O)")):
            acyl_chain = frag
            break
    
    if acyl_chain is None:
        return False, "Failed to isolate acyl chain"
    
    # Check that acyl chain is linear (no rings)
    if Chem.GetSSSR(acyl_chain) > 0:
        return False, "Acyl chain contains rings"
    
    # Identify the carbonyl carbon in the acyl chain fragment
    pattern = Chem.MolFromSmarts("C(=O)[C,H]")
    match = acyl_chain.GetSubstructMatch(pattern)
    if not match:
        return False, "Failed to find carbonyl carbon in acyl chain"
    carbonyl_c_idx = match[0]
    
    # Perform DFS to find the longest linear chain starting from carbonyl carbon
    visited = set()
    max_chain_length = [0]  # Use a list to allow modification within dfs
    
    def dfs(atom_idx, current_length):
        atom = acyl_chain.GetAtomWithIdx(atom_idx)
        visited.add(atom_idx)
        is_terminal = True  # Assume it's a terminal atom until proven otherwise
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in visited:
                continue
            if neighbor.GetAtomicNum() != 6:
                continue  # Only follow carbon atoms
            bond = acyl_chain.GetBondBetweenAtoms(atom_idx, neighbor_idx)
            if bond.IsInRing():
                continue  # Skip bonds in rings
            # Check for branching (degree > 2)
            if neighbor.GetDegree() > 3:
                continue  # Exclude atoms that would introduce branching
            is_terminal = False
            dfs(neighbor_idx, current_length + 1)
        if is_terminal:
            if current_length > max_chain_length[0]:
                max_chain_length[0] = current_length
    
    dfs(carbonyl_c_idx, 1)  # Start length at 1 to include carbonyl carbon
    
    chain_length = max_chain_length[0]
    
    if chain_length > 22:
        return True, f"Acyl chain has length {chain_length}, which is greater than 22"
    else:
        return False, f"Acyl chain has length {chain_length}, which is not greater than 22"