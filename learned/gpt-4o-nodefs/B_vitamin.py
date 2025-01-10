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
        'Thiamine (B1)': Chem.MolFromSmarts('C=C1C([NH3+])=CCN(C)C1=C[CH]S(C)C'),
        'Riboflavin (B2)': Chem.MolFromSmarts('C1(C(C)C[…](OC)c2c(cn3c(nc(=O)[nH]c3=O)n(C)[C@@H]12)C)=O)=O'),
        'Niacin (B3)': Chem.MolFromSmarts('OC(=O)[CH]1CN=CNC1'),
        'Pantothenic acid (B5)': Chem.MolFromSmarts('[C@@H](C(O)=O)NCCC(N)=O'),
        'Pyridoxine (B6)': Chem.MolFromSmarts('C=C1NC([NH3+]CC([C@H]2CO)O)=CC2=C[O]1'),
        'Biotin (B7)': Chem.MolFromSmarts('[C@H]1CS[C@@H](CCCNC(=N)C(=O)[O])N2C(=O)[C@H]12'),
        'Folate (B9)': Chem.MolFromSmarts('N[C@H](CNc1ccc2nc([nH]c3nc(N)[nH]c3=O)c2c1)c1ccc(cc1)[C@@H](C(O)=O)N(C=O)[CH](O)=O'),
        'Cobalamin (B12)': Chem.MolFromSmarts('[C@H]1C[N]2CCC3(C)N4CCC(C=5CNC=6NC=7CNC=8C9NC(C=4C=53)=N1[Co]O(C3N(C=5C=63)[N]7C=9)]=N)C=O')
    }

    # Check for matches
    for vitamin_name, pattern in patterns.items():
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return True, f"Matches {vitamin_name}"
    
    return False, "No B vitamin pattern matched"