"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
"""
Classifies: CHEBI:85021 tetradecanoate ester
"""
from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule contains at least one tetradecanoate (myristoyl) ester group.
    A tetradecanoate ester must have a straight 14-carbon acyl chain (including carbonyl carbon)
    attached via an ester bond with no branching or unsaturation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved ester pattern: O connected to carbon (not hydrogen)
    ester_pattern = Chem.MolFromSmarts('[O;X2]([#6])[C;X3]=[O;X1]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    for match in ester_matches:
        if len(match) < 3:
            continue  # Ensure we have O, C, =O atoms
        
        o_idx, c_idx, _ = match  # O-C=O atoms
        
        # Get the acyl chain carbon (connected to carbonyl, not ester oxygen)
        acyl_carbons = [n for n in mol.GetAtomWithIdx(c_idx).GetNeighbors()
                       if n.GetAtomicNum() == 6 and n.GetIdx() != o_idx]
        
        if not acyl_carbons:
            continue
        
        current = acyl_carbons[0]
        visited = {c_idx, o_idx, current.GetIdx()}
        chain_length = 1  # First carbon after carbonyl
        
        valid_chain = True
        
        while valid_chain and chain_length < 14:
            # Get next carbons (must have exactly 1 non-visited carbon neighbor)
            neighbors = [n for n in current.GetNeighbors()
                        if n.GetAtomicNum() == 6
                        and n.GetIdx() not in visited
                        and mol.GetBondBetweenAtoms(current.GetIdx(), n.GetIdx()).GetBondType() == Chem.BondType.SINGLE]
            
            # Check for branches or end of chain
            if len(neighbors) != 1:
                valid_chain = False
                break
            
            current = neighbors[0]
            visited.add(current.GetIdx())
            chain_length += 1
        
        # Final check: exact length and full traversal
        if valid_chain and chain_length == 14:
            return True, "Contains a straight 14-carbon acyl ester group (myristoyl)"
    
    return False, "No valid myristoyl ester group found"