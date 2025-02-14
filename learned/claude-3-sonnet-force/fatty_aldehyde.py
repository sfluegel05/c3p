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
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for aldehyde functional group (-CHO) at the end of a carbon chain
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde group found"
    
    aldehyde_atom = mol.GetSubstructMatch(aldehyde_pattern)[0]
    aldehyde_carbon = mol.GetAtomWithIdx(aldehyde_atom)
    
    # Check if aldehyde is at the end of the chain
    if aldehyde_carbon.GetDegree() > 2:
        return False, "Aldehyde group not at the end of the chain"
    
    # Check for carbon chain
    carbon_chain = []
    curr_atom = aldehyde_carbon
    while curr_atom.GetAtomicNum() == 6:
        carbon_chain.append(curr_atom.GetIdx())
        neighbors = [nb for nb in curr_atom.GetNeighbors() if nb.GetAtomicNum() == 6]
        if len(neighbors) == 0:
            break
        curr_atom = neighbors[0]
    
    if len(carbon_chain) < 2:
        return False, "No carbon chain found"
    
    return True, "Aldehyde group at the end of a carbon chain"