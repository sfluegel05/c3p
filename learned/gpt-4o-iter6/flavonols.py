"""
Classifies: CHEBI:28802 flavonols
"""
from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is defined as 'Any hydroxyflavone in which the ring hydrogen at
    position 3 of the heterocyclic ring is replaced by a hydroxy group'.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a flavonol, False otherwise
        str: Reason to aid in classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved SMARTS for flavonol skeleton
    # Look for benzopyran-4-one structure with 3-hydroxy group
    flavonol_pattern = Chem.MolFromSmarts('c1cc(O)c2c(c1)[o,O]c(=O)c(c2)O')
    
    # Primary check for flavonol backbone presence
    if not mol.HasSubstructMatch(flavonol_pattern):
        return False, "No flavonol backbone detected"

    # Assess presence of additional hydroxyl groups beyond the core 3-hydroxy
    hydroxy_additions_pattern = Chem.MolFromSmarts('[OH]c1coc2c(c1=O)ccc(c2)O')
    if not mol.HasSubstructMatch(hydroxy_additions_pattern):
        return False, "Essential 3-position hydroxy not identified"

    # Specific check to ensure accurate description of flavonol
    return True, "Structure matches 3-hydroxyflavone backbone consistent with flavonol"

# This logic now reflects potential substitutions while focusing on key structural characteristics of flavonols.