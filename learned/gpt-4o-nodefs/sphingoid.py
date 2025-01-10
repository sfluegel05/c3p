"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    Sphingoids typically have a long aliphatic chain, with characteristic 
    functional groups such as 2-amino-1,3-diol backbones.

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
    
    # Define a SMARTS pattern for sphingoid base - 2-amino-1,3-diol pattern
    # This is a minimal requirement and could vary greatly among different sphingoids
    sphingoid_pattern = Chem.MolFromSmarts("N[C@@H](CO)[C@@H](O)[CX4]")
    
    if mol.HasSubstructMatch(sphingoid_pattern):
        return True, "Contains sphingoid base pattern: 2-amino-1,3-diol"
    
    # Check for long aliphatic chain
    aliphatic_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCC")
    if mol.HasSubstructMatch(aliphatic_chain_pattern):
        return True, "Contains a long aliphatic chain, characteristic of sphingoids"
    
    return None, None  # Placeholder return if classification is ambiguous