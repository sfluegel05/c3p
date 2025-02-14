"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: CHEBI:35748 fatty aldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is an aldehyde formally arising from reduction of the carboxylic acid group
    of its corresponding fatty acid, having a carbonyl group at one end of the carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Preprocess SMILES string to handle backslash characters
    preprocessed_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles.replace('\\', '_')))
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(preprocessed_smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for aldehyde functional group (-CHO) at the end of a carbon chain
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    if not aldehyde_matches:
        return False, "No aldehyde group found"
    
    # Check if the aldehyde group is at the end of a linear carbon chain
    for aldehyde_atom in aldehyde_matches:
        aldehyde_carbon = mol.GetAtomWithIdx(aldehyde_atom)
        
        # Check if aldehyde is at the end of the chain
        if aldehyde_carbon.GetDegree() > 2:
            continue
        
        # Check for linear carbon chain
        carbon_chain = []
        visited = set()
        curr_atom = aldehyde_carbon
        max_chain_length = 30  # Maximum allowed length of the carbon chain
        
        while curr_atom.GetAtomicNum() == 6:
            if curr_atom.GetIdx() in visited:
                break  # Cycle detected, skip this chain
            visited.add(curr_atom.GetIdx())
            carbon_chain.append(curr_atom.GetIdx())
            
            neighbors = [nb for nb in curr_atom.GetNeighbors() if nb.GetAtomicNum() == 6]
            if len(neighbors) != 2:
                break  # Not a linear chain
            
            curr_atom = neighbors[0] if neighbors[0].GetIdx() not in visited else neighbors[1]
            
            if len(carbon_chain) > max_chain_length:
                break  # Chain too long, skip this chain
        
        if len(carbon_chain) >= 2:
            return True, "Aldehyde group at the end of a linear carbon chain"
    
    return False, "No aldehyde group found at the end of a linear carbon chain"