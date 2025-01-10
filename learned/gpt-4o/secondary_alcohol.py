"""
Classifies: CHEBI:35681 secondary alcohol
"""
from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.
    A secondary alcohol has a hydroxyl group (-OH) attached to a saturated carbon
    atom that is bonded to two other carbon atoms, plus allowing configurational specificity.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Attempt to parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Secondary alcohol pattern: Carbon with hydroxyl group having exactly one hydrogen and two carbon neighbors
    secondary_alcohol_pattern = Chem.MolFromSmarts("[CX4H1](C)(C)O")
    
    if mol.HasSubstructMatch(secondary_alcohol_pattern):
        return True, "Contains a hydroxy group attached to a secondary carbon atom"
    
    # If checks fail, it is likely not a secondary alcohol
    return False, "Does not meet the criteria for a secondary alcohol"