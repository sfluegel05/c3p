"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
from rdkit import Chem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    A 17beta-hydroxy steroid must have a hydroxy group at position 17 in beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General pattern for steroid four-ring structure, allowing any stereochemistry
    steroid_pattern = Chem.MolFromSmarts('C1CCC2CCCCC2C1')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Pattern for a hydroxy group at the 17th position in a more general form
    # We use an arbitrary [C@H] for stereochemistry to allow for flexibility in stereochemical description
    hydroxy_17beta_pattern = Chem.MolFromSmarts('[C@H](O)[C@@H]1CCC2C3CCC4CCCC(C3)C4[C@@H]12')
    if not mol.HasSubstructMatch(hydroxy_17beta_pattern):
        return False, "No 17beta-hydroxy group found"

    return True, "Contains 17beta-hydroxy group with steroid backbone configuration"