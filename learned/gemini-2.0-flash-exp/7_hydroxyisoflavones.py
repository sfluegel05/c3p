"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
from rdkit import Chem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    A 7-hydroxyisoflavone has a hydroxyl group at the 7-position of the isoflavone core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 7-hydroxyisoflavone, False otherwise.
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for the 7-hydroxyisoflavone core structure
    # 'X4' is to specify that the hydroxyl group is attached to a carbon with 4 bonds
    core_pattern = Chem.MolFromSmarts('c1cc(O)cc2c1c(=O)oc2')
    if not mol.HasSubstructMatch(core_pattern):
         return False, "Molecule does not have the 7-hydroxyisoflavone core structure"

    return True, "Molecule has the 7-hydroxyisoflavone core with a hydroxyl group at the 7-position"