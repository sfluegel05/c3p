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
    
    # SMARTS pattern for flavonol skeleton with 3-hydroxy group
    flavonol_pattern = Chem.MolFromSmarts('c1cc(O)c2c(c1)[o,O]c(=O)c(c2)O')
    
    # Check for flavonol backbone with 3-hydroxy
    if not mol.HasSubstructMatch(flavonol_pattern):
        return False, "No flavonol backbone detected with 3-hydroxy group"

    # Specific check to ensure accurate description of flavonol
    return True, "Structure matches 3-hydroxyflavone backbone consistent with flavonol"