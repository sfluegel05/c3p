"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
"""
Classifies: CHEBI:47778 pyrimidine deoxyribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside contains a pyrimidine base attached to a deoxyribose sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for pyrimidine base pattern (N1C=CC(=NC1=O))
    pyrimidine_pattern = Chem.MolFromSmarts("[nH]1[c,n][c,n][c,n][c,n]1(=O)")
    if not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No pyrimidine base found"

    # Check for deoxyribose sugar pattern (C1C(O)C(O)C(O1)CO)
    deoxyribose_pattern = Chem.MolFromSmarts("[C@H]1[C@H](O)[C@H](O)[C@@H](CO)O1")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No deoxyribose sugar found"

    # Check for glycosidic bond between pyrimidine base and deoxyribose sugar
    glycosidic_bond_pattern = Chem.MolFromSmarts("[nH]1[c,n][c,n][c,n][c,n]1(=O).[C@H]1[C@H](O)[C@H](O)[C@@H](CO)O1")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond found between pyrimidine base and deoxyribose sugar"

    return True, "Contains pyrimidine base attached to deoxyribose sugar via glycosidic bond"