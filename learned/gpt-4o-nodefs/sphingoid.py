"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    Sphingoids typically feature long aliphatic chains with specific functional groups,
    such as 2-amino-1,3-diol or similar structures.

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

    # Expanded collection of SMARTS patterns for sphingoid recognition
    patterns = [
        # 2-amino-1,3-diol core structure
        Chem.MolFromSmarts("[NX3][C@H](CO)[C@H](O)[CX4]"),
        # Various long-chain configurations with functional groups
        Chem.MolFromSmarts("[CX3,CX4]CC[OX2][C@@H]([NX3])[CX4]"), 
        Chem.MolFromSmarts("[CX3,CX4][CX4][C@H]([NX3])CO"),
        # Double-bond and amino linkage variants
        Chem.MolFromSmarts("C=C[CX4][C@@H](O)[C@@H](N)C"),
        # Broad spectrum pattern for long aliphatic chains with functional groups
        Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][OX2][C@@H](N)[CX4]"),
    ]

    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, f"Matches sphingoid structural pattern: {pattern}"

    # Check for extended aliphatic chains in combination with functional groups
    combined_patterns = [
        Chem.MolFromSmarts("[NX3][C@H]([CX4])[OX2][CX4][CX4][CX4][CX4][CX4][CX4][CX4]"),
        Chem.MolFromSmarts("[C@@H]C=C[CX4][C@@H](O)[C@@H](N)[CX4]"),
    ]

    for pattern in combined_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a combination of long aliphatic chain and functional groups typical for sphingoids."
    
    return False, "Did not match known sphingoid patterns"