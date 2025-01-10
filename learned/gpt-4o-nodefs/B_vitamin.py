"""
Classifies: CHEBI:75769 B vitamin
"""
from rdkit import Chem

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a B vitamin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a B vitamin, False otherwise
        str: Reason for classification, if True indicates type of B vitamin
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define more accurate SMARTS patterns for different B vitamins
    patterns = {
        'Thiamine (B1)': Chem.MolFromSmarts('Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1'),
        'Riboflavin (B2)': Chem.MolFromSmarts('C[C@H]1O[C@@H](CO)n2cnc3c2nc(=O)[nH]c3=O1'),
        'Niacin (B3)': Chem.MolFromSmarts('c1ccccn1C(=O)O'),
        'Pantothenic acid (B5)': Chem.MolFromSmarts('CC(C)(CO)[C@@H](O)C(=O)NCCC(O)=O'),
        'Pyridoxine (B6)': Chem.MolFromSmarts('CC1=C(C)C=C(CO)C(O)=N1'),
        'Biotin (B7)': Chem.MolFromSmarts('N1C(=O)N[C@H]2[C@@H](CCC(=O)O)S[C@@H]12'),
        'Folate (B9)': Chem.MolFromSmarts('Nc1nc2NCC(c3c(=O)[nH]c(nc3)c2c(=O)[nH]1)N'),
        'Cobalamin (B12)': Chem.MolFromSmarts('C[C@]1(CC(N)=O)[C@@H]2C[C@H](OC)[C@H](O)[C@H](O)[C@@H]2O[C@@H]1n3c4cc(C)c(C)cc4c5[n+](c3)c6ccc(nc6)[Co]5')
    }
    
    # Check for matches
    for vitamin_name, pattern in patterns.items():
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return True, f'Matches {vitamin_name}'
    
    return False, "No B vitamin pattern matched"