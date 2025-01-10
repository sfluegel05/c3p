"""
Classifies: CHEBI:35681 secondary alcohol
"""
from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.
    A secondary alcohol possesses a hydroxy group (-OH) attached to a saturated carbon atom
    (sp3 hybridized) which in turn is attached to two other carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alcohol, False otherwise
        str: Reason for classification or failure
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for secondary alcohol
    secondary_alcohol_pattern = Chem.MolFromSmarts("[CX4H](O)[CX4][CX4]")
    
    # Check for pattern in molecule
    if mol.HasSubstructMatch(secondary_alcohol_pattern):
        return True, "Contains a secondary alcohol group"
    else:
        return False, "Does not contain a secondary alcohol group"

# Usage example:
# result, reason = is_secondary_alcohol("CC(O)C")  # Here "CC(O)C" is a classic secondary alcohol structure