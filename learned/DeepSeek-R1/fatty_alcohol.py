"""
Classifies: CHEBI:24026 fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is an aliphatic alcohol with a chain of 3 to over 27 carbons.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for at least one hydroxyl group on an aliphatic carbon
    hydroxyl_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen
            if atom.GetDegree() == 1 and atom.GetTotalNumHs() >= 1:
                carbon = atom.GetNeighbors()[0]
                if carbon.GetAtomicNum() == 6 and not carbon.GetIsAromatic():
                    hydroxyl_found = True
                    break
    if not hydroxyl_found:
        return False, "No hydroxyl group on aliphatic carbon"
    
    # Function to find the longest carbon chain
    def get_longest_chain(mol):
        carbons = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6]
        max_length = 0
        for start in carbons:
            for end in carbons:
                if start == end:
                    continue
                try:
                    path = Chem.GetShortestPath(mol, start, end)
                except:
                    continue
                if len(path) > max_length:
                    max_length = len(path)
        return max_length
    
    chain_length = get_longest_chain(mol)
    if chain_length < 3:
        return False, f"Longest carbon chain is {chain_length} (<3)"
    
    return True, f"Aliphatic alcohol with {chain_length}-carbon chain"