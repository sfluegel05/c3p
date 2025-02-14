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
    (normally adenine, guanine, xanthine, thymine, cytosine, or uracil)
    and either a ribose or deoxyribose sugar as functional parents.

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

    # Define general purine and pyrimidine nucleobase patterns
    purine_base = Chem.MolFromSmarts('c1[nH]c2c([nH]1)ncnc2')  # General purine ring
    pyrimidine_base = Chem.MolFromSmarts('c1c[nH]cnc1')        # General pyrimidine ring

    # Check for presence of nucleobase
    has_nucleobase = False
    if mol.HasSubstructMatch(purine_base):
        has_nucleobase = True
        nucleobase_type = 'purine'
    elif mol.HasSubstructMatch(pyrimidine_base):
        has_nucleobase = True
        nucleobase_type = 'pyrimidine'
    else:
        return False, "No nucleobase (purine or pyrimidine ring) found"

    # Define a general sugar pattern (furanose ring)
    sugar_pattern = Chem.MolFromSmarts('C1OC[C@H](O)[C@@H]1O')  # Five-membered ring with oxygen and hydroxyls

    # Check for presence of sugar ring
    has_sugar = False
    if mol.HasSubstructMatch(sugar_pattern):
        has_sugar = True
    else:
        return False, "No sugar ring found"

    # Define N-glycosidic bond pattern between sugar anomeric carbon and nucleobase nitrogen
    glycosidic_bond_pattern = Chem.MolFromSmarts('[C@H]1([O,N])[C,H][C,H][C,H][O]1[*:1][n,nH][c,c]')  # Anomeric carbon linked to nucleobase nitrogen

    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No N-glycosidic bond between sugar and nucleobase found"

    return True, f"Contains {nucleobase_type} nucleobase attached to sugar via N-glycosidic bond"