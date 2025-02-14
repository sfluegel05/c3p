"""
Classifies: CHEBI:35746 fatty aldehyde
"""
from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is characterized by having an aldehyde group at the end of a carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for the aldehyde pattern: [C]=[O]
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1]=[OX1]")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde group found"
        
    # Check for a likely carbon chain - basic version: ensure at least 4 carbons connected
    carbon_chain_pattern = Chem.MolFromSmarts("C~C~C~C")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Carbon chain is too short"
        
    # Ensure the placement of the aldehyde is at the terminal end of the carbon chain
    # For simplicity, assume a typical linear chain construction with terminal C=O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetTotalDegree() == 1:
            connected_atoms = list(atom.GetNeighbors())
            if len(connected_atoms) == 1 and connected_atoms[0].GetAtomicNum() == 8:
                return True, "Has terminal aldehyde group in long carbon chain"
    
    return False, "Aldehyde group not terminal or improper chain structure"