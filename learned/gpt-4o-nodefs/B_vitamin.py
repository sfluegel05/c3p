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
    
    # Define SMARTS patterns for different B vitamins
    patterns = {
        'Thiamine (B1)': Chem.MolFromSmarts('Cc1nc(N)ncc1C'),
        'Riboflavin (B2)': Chem.MolFromSmarts('C1=CC2=NC(=O)N(C=C2N1)C[C@H](O)[C@H](O)CO'),
        'Niacin (B3)': Chem.MolFromSmarts('n1cccc1C(=O)O'),
        'Pantothenic acid (B5)': Chem.MolFromSmarts('CC(C)(CO)[C@@H](O)C(=O)NCCC(=O)O'),
        'Pyridoxine (B6)': Chem.MolFromSmarts('CC1=NC(C)=C(O1)CO'),
        'Biotin (B7)': Chem.MolFromSmarts('C1C2C(CSC1)NC(=O)NC2'),
        'Folate (B9)': Chem.MolFromSmarts('Nc1nc2nccc(NC3=CC=CC=C3)c2c(=O)[nH]1'),
        'Cobalamin (B12)': Chem.MolFromSmarts('C1C=CC2CC1N1C(C=CC3=[N+]4C(C(CC(N)=O)=C(C3=[N+]6C(CCC(N)=O)=C(C)C4=[N+]5C(CC(N)=O)C(C)C(=C(C)C5=[N+]16)))[N+]72)C1CC(N)=O')
    }
    
    # Check for matches
    for vitamin_name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            return True, f'Matches {vitamin_name}'
    
    return False, "No B vitamin pattern matched"

# Testing with a known B vitamin SMILES
example_smiles = "CN1[C@@H](CNC2=CC=C(C=C2)C(=O)N[C@@H](CCC(O)=O)C(O)=O)CNC2=C1C(=O)NC(N)=N2"
print(is_B_vitamin(example_smiles))