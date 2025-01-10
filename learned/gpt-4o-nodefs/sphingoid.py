"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    Sphingoids generally present a long aliphatic chain and may have characteristic 
    functional groups, such as 2-amino-1,3-diol or variations thereof.

    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Collection of SMARTS patterns for sphingoid recognition
    patterns = [
        Chem.MolFromSmarts("N[C@@H](CO)[C@@H](O)[CX4]"),  # Original 2-amino-1,3-diol backbone
        Chem.MolFromSmarts("N[C@@H](C)[C@@H](O)CO"),      # Variants of amino alcohol structure
        Chem.MolFromSmarts("N[C@H](C)C(O)C=O"),           # Examples with keto presence
        Chem.MolFromSmarts("CC[C@H](O)CO"),               # Part of long-enchainment
        # More SMARTS patterns could be added here based on further research and requirements
    ]

    # Check against each pattern in the collection
    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, f"Matches sphingoid structural pattern: {pattern}"

    # Secondary check: Long aliphatic chain presence
    aliphatic_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC")
    if mol.HasSubstructMatch(aliphatic_chain_pattern):
        return True, "Contains a long aliphatic chain, characteristic of sphingoids"

    return False, "Did not match known sphingoid patterns"