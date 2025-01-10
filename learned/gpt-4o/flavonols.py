"""
Classifies: CHEBI:28802 flavonols
"""
from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    Flavonols are characterized by a 3-hydroxyflavone backbone.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Extended SMARTS pattern to capture core flavonol structure with more variations
    # Flavonol defined by a benzopyranone core with C3 hydroxyl group
    flavonol_pattern = Chem.MolFromSmarts('c1cc(O)c2c(c1)oc(c(=O)c2)-c1ccc(O)cc1')
    
    if not mol.HasSubstructMatch(flavonol_pattern):
        return False, "No flavonol core structure found"

    # Additional criteria for flavonols - more flexible on number/position
    hydroxy_pattern = Chem.MolFromSmarts('[OX2H]')
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)

    if len(hydroxy_matches) < 3:
        return False, f"Insufficient hydroxyl groups, found {len(hydroxy_matches)}"

    # If no issues, confirm it's a flavonol
    return True, "Contains the core structure of a flavonol with sufficient hydroxyl groups"