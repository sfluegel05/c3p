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
        'Pantothenic acid (B5)': Chem.MolFromSmarts('CC(C)(CO)[C@@H](O)C(=O)NCCC(O)=O'),
        'Pyridoxine (B6)': Chem.MolFromSmarts('CC1=NC(C)=C(O1)CO'),
        'Biotin (B7)': Chem.MolFromSmarts('C1C2CS[C@@H](CCCCC(O)=O)[C@@H]2NC(=O)N1'),
        'Folate (B9)': Chem.MolFromSmarts('Nc1nc2n(ccc2c(=O)[nH]1)NCC3=CC=CC=C3'),
        'Cobalamin (B12)': Chem.MolFromSmarts('Oc1nc2c(c(=O)[nH]1)CNC2C')
    }
    
    # Check for matches
    for vitamin_name, pattern in patterns.items():
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return True, f'Matches {vitamin_name}'
    
    return False, "No B vitamin pattern matched"

# Testing with a known B vitamin SMILES
example_smiles = "CN1[C@@H](CNC2=CC=C(C=C2)C(=O)N[C@@H](CCC(O)=O)C(O)=O)CNC2=C1C(=O)NC(N)=N2"
print(is_B_vitamin(example_smiles))