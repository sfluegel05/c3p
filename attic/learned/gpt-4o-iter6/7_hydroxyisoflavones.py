"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    A 7-hydroxyisoflavone possesses the isoflavone skeleton with a hydroxy group specifically at the 7-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 7-hydroxyisoflavone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define isoflavone core pattern (benzopyran-4-one with benzene ring)
    isoflavone_pattern = Chem.MolFromSmarts("Oc1ccc2c(c1)oc(=O)c1ccccc12")
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No isoflavone backbone found"

    # Define 7-hydroxy position SMARTS pattern
    # The position of the hydroxy group is specific, thus defined accordingly
    seven_hydroxy_pattern = Chem.MolFromSmarts("Oc1ccc2c(c1)[o,c](=O)c1ccccc12")
    if not mol.HasSubstructMatch(seven_hydroxy_pattern):
        return False, "No hydroxy group at the 7-position"

    return True, "Has isoflavone skeleton with hydroxy group at 7-position"