"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion has a characteristic CS([O-])(=O)=O group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the alkanesulfonate group pattern: C bonded to S with three oxygens
    sulfonate_pattern = Chem.MolFromSmarts("CS([O-])(=O)=O")
    
    # Check if the molecule matches the sulfonate pattern
    if mol.HasSubstructMatch(sulfonate_pattern):
        return True, "Contains the characteristic alkanesulfonate sulfonate group"
    
    return False, "Does not contain the characteristic alkanesulfonate sulfonate group"