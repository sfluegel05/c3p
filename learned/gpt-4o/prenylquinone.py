"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: Prenylquinone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Define diverse quinone patterns
    quinone_patterns = [
        Chem.MolFromSmarts("C1=CC(=O)C=CC1=O"), # Benzoquinone
        Chem.MolFromSmarts("C1=CC=C(O)C(=O)C=C1"), # Hydroquinone
        Chem.MolFromSmarts("C1=CC=C2C(=C1)C=CC(=O)C2=O"), # Naphthoquinone
        # Additional quinone patterns can be added here
    ]
    
    # Check for the presence of any quinone core structure
    has_quinone = any(mol.HasSubstructMatch(pattern) for pattern in quinone_patterns)
    if not has_quinone:
        return False, "No quinone backbone found"

    # Look for longer prenyl side-chain patterns with repeated isoprenoid units (e.g., C=C-C-C)
    prenyl_pattern_long = Chem.MolFromSmarts("C(=C)CC[C@]+")
    prenyl_matches = mol.GetSubstructMatches(prenyl_pattern_long)
    if len(prenyl_matches) < 2:  # A minimum of 2 matches indicate a decent length side-chain
        return False, "Prenyl side-chain too short"
    
    # Additional checks (e.g., stereochemistry, specific side-chains for subclasses) can be added.

    return True, "Contains quinone backbone with adequate prenyl side-chain"