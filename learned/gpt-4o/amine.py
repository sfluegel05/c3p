"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is derived from ammonia by replacing one, two, or three hydrogen atoms with hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an amine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for primary, secondary, and tertiary amines
    primary_amine_pattern = Chem.MolFromSmarts("[NX3;H2][CX4]")
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3;H1][CX4][CX4]")
    tertiary_amine_pattern = Chem.MolFromSmarts("[NX3][CX4][CX4][CX4]")

    # Check for primary, secondary, or tertiary amine patterns
    if mol.HasSubstructMatch(primary_amine_pattern):
        return True, "Molecule contains a primary amine group"
    elif mol.HasSubstructMatch(secondary_amine_pattern):
        return True, "Molecule contains a secondary amine group"
    elif mol.HasSubstructMatch(tertiary_amine_pattern):
        return True, "Molecule contains a tertiary amine group"

    return False, "No amine group detected"