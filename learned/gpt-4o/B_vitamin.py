"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamins
"""
from rdkit import Chem

def is_B_vitamin(smiles: str):
    """
    Classifies a molecule as a B vitamin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a B vitamin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Define SMARTS patterns for B vitamins
    b_vitamins_patterns = {
        'Vitamin B1 (thiamine)': '[n+]-C-C-S',
        'Vitamin B2 (riboflavin)': 'C1=C(C2=CC(=O)N=C(N=O)N=C2N1)[C@@H](O)C',
        'Vitamin B3 (niacin)': 'nC=O|nC(O)=O',
        'Vitamin B5 (pantothenic acid)': 'CC(C)(CO)C(=O)NCCC(=O)O',
        'Vitamin B6 (pyridoxine)': 'c1c(CO)cnc(C)c1O',
        'Vitamin B7 (biotin)': '[H][C@@]12CS[C@@H](CCCCC(=O)[O-])[C@@]1([H])NC(=O)N2',
        'Vitamin B9 (folic acid)': 'Nc1nc2ncc(CNc3ccc(cc3))nc2c(=O)[nH]1',
        'Vitamin B12 (cobalamin)': 'cobalamin-specific SMILES not easily captured with simple SMARTS'
    }
    
    # Check for substructures
    for name, pattern in b_vitamins_patterns.items():
        substruct = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(substruct):
            return True, f"Matches substructure for {name}"

    return False, "No B vitamin substructure detected"