"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is any oxolane with an oxo-substituent on the tetrahydrofuran ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the tetrahydrofuran ring with an oxo-group pattern
    thf_oxo_pattern = Chem.MolFromSmarts("C1(=O)OC[CH2,CH][CH2,CH]1 |C,C,c,C|")

    # Check for the tetrahydrofuran ring with oxo-group
    if mol.HasSubstructMatch(thf_oxo_pattern):
        return True, "Contains a tetrahydrofuran ring with an oxo group"
    
    return False, "Does not contain a tetrahydrofuran ring with an oxo group"