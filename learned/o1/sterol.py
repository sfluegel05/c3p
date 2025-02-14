"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is any 3β-hydroxy steroid whose skeleton is closely related to cholestan-3-ol,
    allowing for additional carbon atoms in the side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the steroid nucleus (cyclopentanoperhydrophenanthrene)
    steroid_nucleus_smarts = Chem.MolFromSmarts('C1CCC2(CC1)CCC3C2CCC4C3CCCC4')

    if steroid_nucleus_smarts is None:
        return False, "Invalid SMARTS pattern for steroid nucleus"

    # Check for steroid nucleus
    nucleus_match = mol.HasSubstructMatch(steroid_nucleus_smarts)
    if not nucleus_match:
        return False, "No steroid nucleus (cyclopentanoperhydrophenanthrene) found"

    # Define SMARTS pattern for 3β-hydroxyl group attached to the steroid nucleus
    # The 3β-hydroxy group is attached to the third carbon in ring A
    hydroxyl_smarts = Chem.MolFromSmarts('[C@@H]([C;R0])[O;H1]')

    if hydroxyl_smarts is None:
        return False, "Invalid SMARTS pattern for 3β-hydroxyl group"

    # Check for 3β-hydroxyl group
    hydroxyl_match = mol.HasSubstructMatch(hydroxyl_smarts)
    if not hydroxyl_match:
        return False, "No 3β-hydroxyl group found on the steroid backbone"

    return True, "Contains steroid nucleus with 3β-hydroxyl group"