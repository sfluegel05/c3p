"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is an oxolane with an oxo- substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for tetrahydrofuranone: five-membered ring with one oxygen and a carbonyl group
    # The SMARTS pattern accounts for any substituents and a carbonyl group connected
    tetrahydrofuranone_pattern = Chem.MolFromSmarts("O1CC(=O)C(C)CC1")
    
    if mol.HasSubstructMatch(tetrahydrofuranone_pattern):
        return True, "Tetrahydrofuranone structure confirmed"

    return False, "No valid tetrahydrofuranone structure found"