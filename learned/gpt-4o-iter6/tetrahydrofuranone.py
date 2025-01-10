"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is an oxolane having an oxo- substituent at any position
    on the tetrahydrofuran ring.

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

    # SMARTS pattern for a 5-membered lactone ring (tetrahydrofuranone)
    # Modified to capture variations in oxo- substituent location
    lactone_patterns = [
        Chem.MolFromSmarts("O=C1COC[C@H]1"),  # Common oxo-location
        Chem.MolFromSmarts("O=C1COC(C)C1"),   # Dimethyl substitution example
        Chem.MolFromSmarts("O=C1COCC1"),      # Simple tetrahydrofuranone variant
    ]

    # Check for various tetrahydrofuranone structures
    for pattern in lactone_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Tetrahydrofuranone structure confirmed"
    
    return False, "No valid tetrahydrofuranone structure found"