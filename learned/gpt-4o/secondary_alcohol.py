"""
Classifies: CHEBI:35681 secondary alcohol
"""
from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.
    A secondary alcohol has an -OH group attached to an sp3 hybridized carbon, which is bonded to two other carbon atoms.

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

    # Refined SMARTS pattern for a secondary alcohol
    # [C;X4] for sp3 carbon, [OH] for hydroxyl group, connected directly to another [C;X4] indicating secondary carbon
    secondary_alcohol_pattern = Chem.MolFromSmarts("[C;X4](O)[C;X4][C;X4]")
    
    # Check for the presence of the secondary alcohol pattern
    if mol.HasSubstructMatch(secondary_alcohol_pattern):
        return True, "Contains a secondary alcohol functional group"
    
    return False, "No secondary alcohol functional group found"