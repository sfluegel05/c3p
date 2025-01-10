"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
from rdkit import Chem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string,
    requiring a steroid backbone and a stereospecific 11beta-hydroxy group.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more generalized steroid backbone pattern
    # Pattern attempts to capture the [6-6-6-5] steroid nucleus with some flexibility
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3CCCC4")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "Steroid backbone not found"

    # Check for 11beta-hydroxy group
    # Assumes the presence of an -OH group in a beta configuration at a possible 11 position
    # The configuration is defined via stereochemistry descriptors post attachment
    hydroxy_11beta_pattern = Chem.MolFromSmarts("[C@@H]1(O)CCC2C1CCC3C2CCC=C3")
    if not mol.HasSubstructMatch(hydroxy_11beta_pattern):
        return False, "No 11beta-hydroxy group with correct configuration found"

    return True, "Contains steroid backbone with an 11beta-hydroxy group of the expected configuration"