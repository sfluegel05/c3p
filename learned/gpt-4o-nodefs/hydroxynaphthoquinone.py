"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    A hydroxynaphthoquinone contains a naphthoquinone core with one or more hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a hydroxynaphthoquinone, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for naphthoquinone core (two keto groups on a naphthalene system)
    naphthoquinone_pattern = Chem.MolFromSmarts("O=C1C=CC2=CC=CC=C2C1=O")  # Templated structure
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")  # Hydroxyl groups

    # Check for naphthoquinone core
    if not mol.HasSubstructMatch(naphthoquinone_pattern):
        return False, "No recognized naphthoquinone core found"
    
    # Check for at least one hydroxyl group
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) == 0:
        return False, "No hydroxyl group found"
    
    return True, "Contains naphthoquinone core with at least one hydroxyl group"