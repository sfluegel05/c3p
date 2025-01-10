"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: sterol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol
    (additional carbon atoms may be present in the side chain).
    
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
    
    # Define steroid backbone pattern (cyclopentanoperhydrophenanthrene nucleus)
    steroid_pattern = Chem.MolFromSmarts('''
        [#6]1~[#6]~[#6]2~[#6](~[#6]1)~[#6]3~[#6](~[#6]~[#6]~[#6]~3)~[#6]~2
    ''')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Does not contain steroid backbone"
    
    # Define 3-hydroxy group on ring A
    hydroxy_pattern = Chem.MolFromSmarts('[#6]-1(-[#8H])-[#6]=[#6]-[#6]-[#6]-1')
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "Does not have 3-hydroxy group at position 3"
    
    return True, "Contains steroid backbone with 3-hydroxy group characteristic of sterols"

__metadata__ = {   'chemical_class': {   'name': 'sterol',
                              'definition': 'Any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol (additional carbon atoms may be present in the side chain).'}}