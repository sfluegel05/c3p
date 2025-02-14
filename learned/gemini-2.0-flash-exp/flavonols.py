"""
Classifies: CHEBI:28802 flavonols
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is a flavone with a hydroxyl group at the 3-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavone core pattern using SMARTS, allowing for substitutions
    # The core is c1ccccc1-c2oc(c3ccccc3)c(=O)c2 but with substitutions allowed on the rings (not hydrogens)
    flavone_core = Chem.MolFromSmarts('c1ccc2oc(-c3ccccc3)c(=O)c2c1')
    if not mol.HasSubstructMatch(flavone_core):
        return False, "Not a flavone core structure"

    # Check for the hydroxyl group at position 3 of the pyran ring (C ring).
    # This is the carbon directly attached to the carbonyl group in the flavone core
    hydroxyl_at_3 = Chem.MolFromSmarts('c1ccccc1-c2oc(-[C;H1]([O;H1])c(=O)c3ccccc3)c2')
    if not mol.HasSubstructMatch(hydroxyl_at_3):
      return False, "No hydroxyl group at position 3 of the flavone core"

    # If both the flavone core and 3-OH are found, it's a flavonol
    return True, "Flavonol: Flavone core with hydroxyl at position 3"