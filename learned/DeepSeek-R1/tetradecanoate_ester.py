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
    
    # Find all ester groups (O-C=O), explicit valence matching
    ester_pattern = Chem.MolFromSmarts('[O;X2]-[C;X3]=[O;X1]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    for match in ester_matches:
        # Match returns (O, C, =O) indices, we need O and C positions
        if len(match) < 2:
            continue
        o_idx, c_idx = match[0], match[1]
        
        # Get the carbon adjacent to carbonyl (acyl chain direction)
        neighbors = [n for n in mol.GetAtomWithIdx(c_idx).GetNeighbors() 
                    if n.GetIdx() != o_idx and n.GetAtomicNum() == 6]
        if not neighbors:
            continue
        current_atom = neighbors[0]
        
        # Traverse the carbon chain
        chain_length = 1  # Already counting the first carbon after carbonyl
        visited = {c_idx, o_idx, current_atom.GetIdx()}
        
        while True:
            # Follow only single bonds to carbons
            next_carbons = [n for n in current_atom.GetNeighbors()
                           if n.GetAtomicNum() == 6
                           and n.GetIdx() not in visited
                           and mol.GetBondBetweenAtoms(current_atom.GetIdx(), n.GetIdx()).GetBondType() == Chem.BondType.SINGLE]
            
            if len(next_carbons) != 1:
                break  # Branch or end of chain
            current_atom = next_carbons[0]
            visited.add(current_atom.GetIdx())
            chain_length += 1
        
        # Check chain length (13 carbons after carbonyl = total 14)
        if chain_length == 13:
            return True, "Contains a myristoyl (14-carbon) ester group"
    
    return False, "No 14-carbon ester group found"