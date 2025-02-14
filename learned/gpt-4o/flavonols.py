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

    # Extended SMARTS pattern to capture core flavonol structure
    # This pattern matches the 3-hydroxyflavone structure with potential modifications
    flavonol_pattern = Chem.MolFromSmarts('Oc1cc(O)c2c(c1)oc(-c1ccccc1)c(=O)c2')
    
    if not mol.HasSubstructMatch(flavonol_pattern):
        return False, "No flavonol core structure found"

    # Check for variations such as different potential hydroxy positions
    hydroxy_pattern = Chem.MolFromSmarts('[OX2H]')
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)

    if len(hydroxy_matches) < 3:
        return False, f"Insufficient hydroxyl groups, found {len(hydroxy_matches)}"

    return True, "Contains the core structure of a flavonol with sufficient hydroxyl groups"