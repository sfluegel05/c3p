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
    steroid_nucleus_smarts = Chem.MolFromSmarts('''
    [#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1
    ~[#6]2~[#6]~[#6]~[#6]~[#6]~1~[#6]3~[#6]~[#6]~[#6]~[#6]~[#6]3~[#6]~2
    ''')

    # Check for steroid nucleus
    nucleus_match = mol.HasSubstructMatch(steroid_nucleus_smarts)
    if not nucleus_match:
        return False, "No steroid nucleus (cyclopentanoperhydrophenanthrene) found"

    # Define SMARTS pattern for 3β-hydroxyl group attached to the steroid nucleus
    hydroxyl_smarts = Chem.MolFromSmarts('''
    [#6]-1(-[O;H1])-[#6]=[#6]-[#6]-[#6]-[#6]-1
    ''')

    # Check for 3β-hydroxyl group
    hydroxyl_match = mol.HasSubstructMatch(hydroxyl_smarts)
    if not hydroxyl_match:
        return False, "No 3β-hydroxyl group found on the steroid backbone"

    return True, "Contains steroid nucleus with 3β-hydroxyl group"