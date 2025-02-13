"""
Classifies: CHEBI:35681 secondary alcohol
"""
from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.
    A secondary alcohol has an -OH group attached to a saturated carbon, which is bonded to two other carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for secondary alcohol
    # [#6] represents any carbon, [CX4] ensures the carbon is sp3 hybridized (saturated)
    secondary_alcohol_pattern = Chem.MolFromSmarts("[#6][CX4!H0][CH](O)[CX4!H0]")

    # Check for the presence of the secondary alcohol pattern
    if mol.HasSubstructMatch(secondary_alcohol_pattern):
        return True, "Contains a secondary alcohol functional group"
    
    return False, "No secondary alcohol functional group found"