"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
"""
Classifies: 7-hydroxyisoflavones
"""

from rdkit import Chem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    A 7-hydroxyisoflavone is a hydroxyisoflavone compound having a hydroxy group at the 7-position.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a 7-hydroxyisoflavone, False otherwise.
        str: Reason for classification.
    """

    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the isoflavone core with atom mapping
    isoflavone_core_smarts = """
        [#6]-1(:[#6]:[#6]:[#6]:[#6]:[#6]:1)-[#6]-2=[#8]-[#6]3=[#6]([#6]:[#6]:[#6]:[#6]:3)-[#8]=[#6]-2
    """
    isoflavone_core_smarts = isoflavone_core_smarts.replace('\n', '')
    isoflavone_core = Chem.MolFromSmarts(isoflavone_core_smarts)
    if isoflavone_core is None:
        return False, "Invalid SMARTS pattern for isoflavone core"

    # Define the SMARTS pattern for the 7-hydroxy group
    # Assuming position 7 is the carbon in ring A adjacent to the oxygen in the pyran ring
    hydroxy_at_7_smarts = """
        [#6]-1(-[#8]-[#1]):[#6]:[#6]:[#6]:[#6]:[#6]:1
    """
    hydroxy_at_7_smarts = hydroxy_at_7_smarts.replace('\n', '')
    hydroxy_at_7 = Chem.MolFromSmarts(hydroxy_at_7_smarts)
    if hydroxy_at_7 is None:
        return False, "Invalid SMARTS pattern for hydroxy at position 7"

    # Check for isoflavone core
    if not mol.HasSubstructMatch(isoflavone_core):
        return False, "Molecule does not contain the isoflavone core"

    # Check for hydroxy group at position 7
    if not mol.HasSubstructMatch(hydroxy_at_7):
        return False, "Molecule does not have a hydroxy group at position 7"

    return True, "Molecule is a 7-hydroxyisoflavone"