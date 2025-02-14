"""
Classifies: CHEBI:50753 isoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is a 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Define the correct 1-benzopyran core
    # The core pattern:
    core_pattern = Chem.MolFromSmarts("c1ccccc2[O;X2][CH2][C]([*:1])[c]21")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Not a 1-benzopyran core"

    # 2. Define a general aryl group pattern
    aryl_pattern = Chem.MolFromSmarts("c1ccccc1")

    # 3. Create a combined SMARTS pattern that enforces the connection between the benzopyran and the aryl
    # and that this aryl is directly bound to the 3 position of the benzopyran,
    # allowing for a single or double bond
    combined_pattern = Chem.MolFromSmarts("c1ccccc2[O;X2][CH2][C](~c3ccccc3)[c]21")
    
    if not mol.HasSubstructMatch(combined_pattern):
        return False, "No aryl group at position 3"

    return True, "1-Benzopyran with an aryl substituent at position 3"