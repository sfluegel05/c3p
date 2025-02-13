"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol is defined as an alcohol derived from an aliphatic compound,
    with at least one hydroxyl group attached to a non-aromatic carbon (sp3 hybridized).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for hydroxyl group attached to non-aromatic carbon
    hydroxyl_pattern = Chem.MolFromSmarts('[CX4;!H0][OH]')  # Ensure it's attached to an sp3 carbon that is not part of an aromatic system

    # Check for presence of the pattern
    if mol.HasSubstructMatch(hydroxyl_pattern):
        return True, "Molecule contains a hydroxyl group attached to an aliphatic (sp3) carbon, indicative of an aliphatic alcohol"
    else:
        return False, "No suitable hydroxyl group pattern for an aliphatic alcohol found"