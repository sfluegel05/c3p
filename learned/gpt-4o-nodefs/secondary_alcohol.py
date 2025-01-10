"""
Classifies: CHEBI:35681 secondary alcohol
"""
from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.
    A secondary alcohol has a hydroxyl group (-OH) attached to a carbon atom, which is also connected to two other carbon atoms.

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

    # Define SMARTS pattern for a secondary alcohol (O must be attached to a C with 2 other carbon neighbors)
    secondary_alcohol_pattern = Chem.MolFromSmarts("[C;H1,H2][O][C;H2]")
    if mol.HasSubstructMatch(secondary_alcohol_pattern):
        return True, "Contains a secondary alcohol group (R2CHOH structure)"
    
    return False, "Does not contain a secondary alcohol group"