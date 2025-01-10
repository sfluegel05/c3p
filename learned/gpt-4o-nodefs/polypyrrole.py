"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    Since the definition is ambiguous, this analysis extends to multiple pyrrole rings
    and potential polymeric structures containing pyrroles.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as polypyrrole, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Consider pyrrole and possible linked multiple pyrrole structures
    pyrrole_pattern = Chem.MolFromSmarts('n1cccc1')
    
    # If multiple pyrroles might be required in sequence, look for linked dimers or more
    polypyrrole_pattern = Chem.MolFromSmarts('n1cccc1-n2cccc2')
    
    # Check if any pattern relevant to polypyrrole structures matches
    if mol.HasSubstructMatch(polypyrrole_pattern):
        return True, "Contains multiple linked pyrrole rings pattern"

    # Fallback on single pyrrole detection if polymeric form not specific
    if mol.HasSubstructMatch(pyrrole_pattern):
        return True, "Contains pyrrole ring, but not polymeric"

    return False, "No polypyrrole-defining structure found or definition is too vague"