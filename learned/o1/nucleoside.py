"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: CHEBI:33838 nucleoside
"""
from rdkit import Chem

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside is an N-glycosyl compound that has both a nucleobase
    (adenine, guanine, xanthine, thymine, cytosine, or uracil) and either
    a ribose or deoxyribose sugar as functional parents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define nucleobases
    nucleobases = {
        'adenine': Chem.MolFromSmarts('c1ncnc2ncnc12'),
        'guanine': Chem.MolFromSmarts('c1nc2c(n1)[nH]c(=O)[nH]2'),
        'xanthine': Chem.MolFromSmarts('O=C1NC(=O)NC=N1'),
        'thymine': Chem.MolFromSmarts('O=C1NC(=O)C=CN1C'),
        'cytosine': Chem.MolFromSmarts('O=C1NC=CN=C1N'),
        'uracil': Chem.MolFromSmarts('O=C1NC=CC(=O)[nH]1')
    }

    # Check for presence of nucleobase
    has_nucleobase = False
    for name, pattern in nucleobases.items():
        if mol.HasSubstructMatch(pattern):
            has_nucleobase = True
            nucleobase_name = name
            break
    if not has_nucleobase:
        return False, "No nucleobase (adenine, guanine, xanthine, thymine, cytosine, or uracil) found"

    # Define sugar patterns (ribose and deoxyribose)
    # Ribose: a five-membered ring with oxygen and four hydroxyl groups
    ribose = Chem.MolFromSmarts('C1(CO)[C@@H](O)[C@H](O)[C@@H](O)O1')
    # Deoxyribose: similar but missing one hydroxyl group
    deoxyribose = Chem.MolFromSmarts('C1(CO)[C@@H](O)[C@H](O)[C@@H](O)O1')

    # Check for presence of sugar
    has_sugar = False
    if mol.HasSubstructMatch(ribose):
        has_sugar = True
        sugar_name = 'ribose'
    elif mol.HasSubstructMatch(deoxyribose):
        has_sugar = True
        sugar_name = 'deoxyribose'
    else:
        return False, "No ribose or deoxyribose sugar found"

    # Check for N-glycosidic bond between nucleobase and sugar
    # The bond is between an anomeric carbon of sugar and a nitrogen atom of nucleobase
    # Define SMARTS for N-glycosidic bond pattern
    glycosidic_bond = Chem.MolFromSmarts('[nX2]-[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O')
    if not mol.HasSubstructMatch(glycosidic_bond):
        return False, "No N-glycosidic bond between nucleobase and sugar found"

    return True, f"Contains {nucleobase_name} nucleobase attached to {sugar_name} via N-glycosidic bond"