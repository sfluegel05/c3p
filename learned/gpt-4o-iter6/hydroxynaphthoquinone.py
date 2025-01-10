"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    A hydroxynaphthoquinone is a naphthoquinone moiety with at least one hydroxy group substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxynaphthoquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns
    naphthoquinone_pattern = Chem.MolFromSmarts('c1cc2ccccc2c(=O)c(=O)c1')
    hydroxy_pattern = Chem.MolFromSmarts('[OH]')
    
    # Check for naphthoquinone structure
    if not mol.HasSubstructMatch(naphthoquinone_pattern):
        return False, "No naphthoquinone core found"
    
    # Check for at least one hydroxy group
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No hydroxy group found"
    
    # Verify hydroxyl group is attached to naphthoquinone
    # Requires hydroxy group directly connected to naphthoquinone aromatic carbons
    for match in mol.GetSubstructMatches(naphthoquinone_pattern):
        atoms = set(match)
        hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
        if any(h[0] in atoms for h in hydroxy_matches):
            return True, "Contains naphthoquinone core with hydroxy group substitution"
    
    return False, "Hydroxy group not attached to naphthoquinone"