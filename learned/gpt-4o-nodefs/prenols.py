"""
Classifies: CHEBI:26244 prenols
"""
from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    Prenols generally consist of multiple isoprene units and contain an alcohol group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of alcohol group (-OH)
    alcohol_pattern = Chem.MolFromSmarts('[OX2H]')
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "Missing alcohol group"

    # Enhanced check for isoprene units; we broadened the pattern to capture alternative connections
    # Isoprene units: (C=C-C-C) or (C-C=C-C) or (C=C-C) or (C-C=C)
    isoprene_patterns = [
        Chem.MolFromSmarts('C=C-C-C'),
        Chem.MolFromSmarts('C-C=C-C'),
        Chem.MolFromSmarts('C=C-C'),
        Chem.MolFromSmarts('C-C=C')
    ]

    isoprene_count = 0
    for pattern in isoprene_patterns:
        matches = mol.GetSubstructMatches(pattern)
        isoprene_count += len(matches)

    if isoprene_count < 1:
        return False, f"Insufficient isoprene-like units (found {isoprene_count})"

    return True, "Contains alcohol group and isoprene units, typical of prenols"