"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone, or lactone, is an oxolane with an oxo- substituent.

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
    
    # SMARTS pattern for a 5-membered lactone ring with an oxo- substituent
    # This pattern represents a generic lactone: [O]C1CCCC1=O, also capturing tetrahydrofuranone structures
    lactone_pattern = Chem.MolFromSmarts("O=C1OC(C)C1")
    
    # Check for tetrahydrofuranone structure
    if mol.HasSubstructMatch(lactone_pattern):
        return True, "Tetrahydrofuranone structure confirmed"
    
    return False, "No valid tetrahydrofuranone structure found"