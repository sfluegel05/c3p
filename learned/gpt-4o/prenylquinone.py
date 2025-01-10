"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: Prenylquinone
"""
from rdkit import Chem

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone is a quinone substituted by a polyprenyl-derived side-chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define quinone patterns
    quinone_patterns = [
        Chem.MolFromSmarts("C1=CC(=O)C=CC1=O"),            # Benzoquinone
        Chem.MolFromSmarts("C1=CC=C(C=O)C=C1O"),           # Hydroquinone
        Chem.MolFromSmarts("C1=CC=C2C(=C1)C=CC(=O)C2=O"),  # Naphthoquinone
        # Additional specific patterns can be added
    ]

    # Check for the presence of any quinone core structure
    has_quinone = any(mol.HasSubstructMatch(pattern) for pattern in quinone_patterns)
    if not has_quinone:
        return False, "No quinone backbone found"

    # Define prenyl side-chain patterns (e.g., repeated isoprenoid units)
    isoprene_unit = Chem.MolFromSmarts("C(=C)CC")
    prenyl_matches = mol.GetSubstructMatches(isoprene_unit)

    # We expect multiple isoprene units for a significant prenyl side-chain
    if len(prenyl_matches) < 2:
        return False, "Prenyl side-chain too short"

    return True, "Contains quinone backbone with adequate prenyl side-chain"