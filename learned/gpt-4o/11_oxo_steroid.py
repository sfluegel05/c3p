"""
Classifies: CHEBI:47787 11-oxo steroid
"""
from rdkit import Chem

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    An 11-oxo steroid has a ketone group at the 11th carbon position inherited from a steroid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a flexible 11-oxo pattern for steroids.
    # Including aromatic steroids and variations on the traditional A/B/C/D rings.
    # Steroid backbone with oxo on position C11 could be more structurally flexible:
    oxo_pattern = Chem.MolFromSmarts(
        "[#6]1(-[#6]-[#8]=O)(-[#6]-[#6]-[#6])-[#6]-[#6]-[#6]2-[#6]-[#8]=O"
    )
    
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "11-oxo group not found in standard positions"
    
    return True, "Contains 11-oxo group in the correct steroid configuration"