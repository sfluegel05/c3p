"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    Amines are compounds containing a nitrogen atom bonded to hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Smart patterns for primary, secondary, and tertiary amines
    primary_amine_pattern = Chem.MolFromSmarts("[NX3H2]")
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3H1]([#6])")
    tertiary_amine_pattern = Chem.MolFromSmarts("[NX3]([#6])([#6])")

    # Check for presence of primary, secondary, or tertiary amine
    if mol.HasSubstructMatch(primary_amine_pattern):
        return True, "Contains a primary amine group"
    elif mol.HasSubstructMatch(secondary_amine_pattern):
        return True, "Contains a secondary amine group"
    elif mol.HasSubstructMatch(tertiary_amine_pattern):
        return True, "Contains a tertiary amine group"
    
    return False, "No amine group found in the molecule"