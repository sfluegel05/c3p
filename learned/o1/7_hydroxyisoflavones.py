"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
"""
Classifies: 7-hydroxyisoflavones
"""

from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Define the SMARTS pattern for the isoflavone core with a hydroxy group at position 7
    # This pattern represents the isoflavone core and a hydroxy group attached to position 7
    pattern_smarts = '[OH]c1ccc2c(c1)cc(=O)oc2-c1ccccc1'  # Hydroxy at position 7

    # Create the pattern molecule from SMARTS
    pattern = Chem.MolFromSmarts(pattern_smarts)
    if pattern is None:
        return False, "Invalid SMARTS pattern for 7-hydroxyisoflavone core"

    # Check if the molecule matches the 7-hydroxyisoflavone pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Molecule contains the 7-hydroxyisoflavone core with hydroxy at position 7"
    else:
        return False, "Molecule does not contain the 7-hydroxyisoflavone core with hydroxy at position 7"