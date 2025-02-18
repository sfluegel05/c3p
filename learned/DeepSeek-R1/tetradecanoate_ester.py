"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
"""
Classifies: CHEBI:85021 tetradecanoate ester
"""
from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is an ester derived from tetradecanoic acid (myristic acid),
    meaning it contains a myristoyl group (14-carbon chain) attached via an ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all ester groups (O-C=O)
    ester_pattern = Chem.MolFromSmarts('[O]-[C]=O')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    for match in ester_matches:
        o_idx, c_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(c_idx)
        
        # Get the carbon adjacent to the carbonyl (part of the acyl chain)
        acyl_start = [n for n in carbonyl_atom.GetNeighbors() if n.GetIdx() != o_idx and n.GetAtomicNum() == 6]
        if not acyl_start:
            continue
        acyl_start = acyl_start[0]
        
        # Traverse the acyl chain to check length
        current_atom = acyl_start
        chain_length = 1
        visited = {c_idx, o_idx, current_atom.GetIdx()}
        
        while True:
            next_carbons = [n for n in current_atom.GetNeighbors() 
                            if n.GetAtomicNum() == 6 
                            and n.GetIdx() not in visited 
                            and n.GetBond(current_atom).GetBondType() == Chem.BondType.SINGLE]
            
            if len(next_carbons) != 1:
                break  # Branch or end of chain
            current_atom = next_carbons[0]
            visited.add(current_atom.GetIdx())
            chain_length += 1
        
        # Check if chain ends with a methyl group (degree 1) and length is 13 (14 carbons total)
        if current_atom.GetDegree() == 1 and chain_length == 13:
            return True, "Contains a myristoyl (tetradecanoyl) ester group"
    
    return False, "No myristoyl ester group found"