"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.

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

    # Enhanced steroid backbone check using a more flexible pattern
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3CCC(C4)C")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "Basic steroid backbone not found"

    # Identify 11th carbon in steroid backbone to check for beta-OH
    core_match = mol.GetSubstructMatch(steroid_backbone_pattern)
    if not core_match:
        return False, "Core steroid scaffold not matching"
    
    # Check around 11th carbon for -OH group with beta specific configuration
    # The index of the 11th carbon might vary based on match; here using:
    hydroxy_11beta = Chem.MolFromSmarts("C[C@@H](O)")
    has_hydroxy = mol.HasSubstructMatch(hydroxy_11beta)
    
    if not has_hydroxy:
        return False, "No 11beta-hydroxy group with correct configuration found"

    return True, "Contains steroid backbone with an 11beta-hydroxy group of the expected configuration"