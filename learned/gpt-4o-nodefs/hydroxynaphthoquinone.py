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
    
    # Define the SMARTS pattern for naphthoquinone core with flexible positions of carbonyls
    naphthoquinone_patterns = [
        Chem.MolFromSmarts("O=C1C=CC2=CC=CC=C2C1=O"),  # Basic pattern, starting from keto groups
        Chem.MolFromSmarts("C1=CC2=CC(C=CC2C=C1)=O"),  # Adjusted to allow aroma
        Chem.MolFromSmarts("C1=C(O)C=CC2=CC=CC=C2C1=O"), # Tautomer considerations
    ]
    
    # Check if naphthoquinone core is present
    core_match_found = any(mol.HasSubstructMatch(pattern) for pattern in naphthoquinone_patterns)
    if not core_match_found:
        return False, "No recognized naphthoquinone core found"

    # Define the SMARTS pattern for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    
    # Check for at least one hydroxyl group
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) == 0:
        return False, "No hydroxyl group found"

    return True, "Contains naphthoquinone core with at least one hydroxyl group"