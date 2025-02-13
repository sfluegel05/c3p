"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    A 3beta-hydroxy steroid is characterized by a steroid structure with a hydroxyl
    group at the 3rd carbon in the beta configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for steroid backbone pattern (basic tetracyclic structure of steroids)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3C1CCC4C3CCC2(C4)")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Look for 3beta-hydroxy group (3rd carbon with beta stereo OH group)
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H](O)C1CCC2C3C1CCC4C3CCC2(C4)")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3beta-hydroxy group found"

    return True, "Contains steroid backbone with 3beta-hydroxy group"