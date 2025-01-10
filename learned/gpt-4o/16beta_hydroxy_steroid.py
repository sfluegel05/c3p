"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
from rdkit import Chem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    A 16beta-hydroxy steroid is defined as a steroid with an -OH group at position 16 in beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General steroid backbone pattern: four fused rings (commonly denoted as ABCD ring system)
    # C1-C2-C3-C4 core represents simplified steroid structure
    # A closer match reduces specificity, consider known steroid patterns
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CC2CC3C(C1)C2CC[C@@H]3")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No general steroid backbone found"
    
    # Specific substructure matching 16-beta-hydroxy group 
    # Need to specify stereochemistry using RDKit SMARTS stereo markers
    beta_oh_16_pattern = Chem.MolFromSmarts("C[C@H](O)[C@]C")
    if mol.HasSubstructMatch(beta_oh_16_pattern):
        return True, "Contains 16beta-hydroxy group with correct steroid backbone"
    
    return False, "Does not match 16beta-hydroxy steroid criteria"