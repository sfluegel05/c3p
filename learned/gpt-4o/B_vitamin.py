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
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for B vitamins
    b_vitamins_patterns = {
        'Vitamin B1 (thiamine)': 'Cn[c+](cSC)CCOP([O-])([O-])=O',  # Extended SMILES for specificity
        'Vitamin B2 (riboflavin)': 'CNC1=NC2=C(NC1=O)C(=C(N=C2O)N)c3cnc4c(c34)O',
        'Vitamin B3 (niacin and derivatives)': 'c1cc(cnc1)C(=O)O',
        'Vitamin B5 (pantothenic acid and derivatives)': 'C(C(C(C(=O)NCCCO)C(=O)O)O)C',
        'Vitamin B6 (pyridoxine and derivatives)': 'c1cc(COP([O-])([O-])=O)c(O)cn1C',
        'Vitamin B7 (biotin)': '[H][C@@]12CSCC(=O)N1[C@H](CCC(=O)[O-])C2',
        'Vitamin B9 (folic acid and derivatives)': 'Nc1nc2ncc(NC(=O)C3=CC=C(C=C3)C(=O)NC[C@@H](CCC([O-])=O)C([O-])=O)c2[nH]1',
        'Vitamin B12 (cobalamin-related structures)': '[Co]N5C(=C1C=C4[C@@H]2[C@H]3N1C=CC25)C3=CC(=C4)C'
    }
    
    # Check for substructures
    for name, pattern in b_vitamins_patterns.items():
        substruct = Chem.MolFromSmarts(pattern)
        if not substruct:
            continue
        if mol.HasSubstructMatch(substruct):
            return True, f"Matches substructure for {name}"

    return False, "No B vitamin substructure detected"