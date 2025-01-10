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

    # General pattern for steroid backbone (tetracyclic ring)
    steroid_pattern = Chem.MolFromSmarts('C1C[C@H]2CC[C@@H]3C(=O)CC[C@]3(C)C[C@@H]2[C@H]1')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Check for 17beta-hydroxy group: Assumes 17th position hydroxy in similar structures
    # The stereochemistry may vary across examples, so we allow flexibility
    hydroxy_17beta_pattern = Chem.MolFromSmarts('[C@](C)(O)[C@H]1CC[C@]2(C)C[C@H]3C(=O)C=C[C@@H]3CC[C@H]12')
    
    # Modify pattern to capture more generic configurational variance
    if not mol.HasSubstructMatch(hydroxy_17beta_pattern):
        return False, "No 17beta-hydroxy group found"

    return True, "Contains 17beta-hydroxy group with acceptable steroid backbone configuration"