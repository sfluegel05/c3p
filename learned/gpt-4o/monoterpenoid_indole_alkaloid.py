"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Define the SMARTS pattern for an indole ring
    indole_pattern = Chem.MolFromSmarts('c1ccc2nc3ccccc3c2c1')
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole ring found"

    # Define a broad pattern for a terpenoid portion via isoprene units (C5 isoprene structure)
    terpenoid_pattern = Chem.MolFromSmarts('C=C(C)C')
    if not mol.HasSubstructMatch(terpenoid_pattern):
        return False, "No terpenoid features found (missing isoprene units)"
    
    # Additional checks involving molecular weight or specific atom counts could be added if needed

    # If passes all key structure tests, classify as a monoterpenoid indole alkaloid
    return True, "Contains both indole ring and terpenoid features (isoprene units)"