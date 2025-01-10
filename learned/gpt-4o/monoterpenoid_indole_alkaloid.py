"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
from rdkit import Chem

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    A monoterpenoid indole alkaloid is characterized by a structure containing an indole ring
    along with terpenoid-like features typically originating from a diisoprenoid unit.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expanded SMARTS pattern for an indole ring to handle various substructures
    indole_pattern = Chem.MolFromSmarts('[nH]1c2ccccc2c3ccccc13')
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole ring found"

    # Expanded SMARTS pattern for terpenoid-like features
    # A general example pattern for the types commonly found in biological terpenes
    terpenoid_pattern = Chem.MolFromSmarts('C(C)(C)C=C')  # More complex than a simple isoprene
    if not mol.HasSubstructMatch(terpenoid_pattern):
        return False, "No terpenoid features found (missing complex isoprene-related structures)"
    
    # If passes all key structure tests, classify as a monoterpenoid indole alkaloid
    return True, "Contains both indole ring and terpenoid features."