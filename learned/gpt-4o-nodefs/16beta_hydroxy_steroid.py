"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.

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

    # Define steroid backbone: a typical perhydrocyclopentanophenanthrene system
    steroid_backbone = Chem.MolFromSmarts('C1CC2C3C=C4CCC3C2C1CCC4')
    if not mol.HasSubstructMatch(steroid_backbone):
        return False, "No steroid backbone found"

    # Identify the hydroxy group at the 16beta position using SMARTS
    # This is a simplification, assuming correct stereochemistry is determined externally
    # for just the presence of a hydroxyl. In practice, this would also verify stereochemistry.
    hydroxyl_16beta_pattern = Chem.MolFromSmarts('C3C[C@@H](O)C2CCC1C[C@@H]1C2')
    
    # Ensure the molecule has the 16beta hydroxyl pattern
    if not mol.HasSubstructMatch(hydroxyl_16beta_pattern):
        return False, "No 16beta-hydroxy group found"

    return True, "Contains steroid backbone with 16beta-hydroxy group"