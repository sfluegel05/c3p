"""
Classifies: CHEBI:26255 prenylquinone
"""
from rdkit import Chem

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone is a quinone substituted by a polyprenyl-derived side chain.
    
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
    
    # Define SMARTS pattern for quinone core
    quinone_pattern = Chem.MolFromSmarts("C1=CC(=O)[C@H]=CC1=O")  # Simplified example
    
    if not mol.HasSubstructMatch(quinone_pattern):
        return False, "No quinone structure found"
        
    # Define SMARTS pattern for prenyl chain (repeated isoprene units)
    # Simplified version: it looks for alkenyl chains indicating polyprenyl-like structure
    prenyl_chain_pattern = Chem.MolFromSmarts("C=C-C=C")  # Basic example
    
    # Check for polyprenyl side chain
    if not mol.HasSubstructMatch(prenyl_chain_pattern):
        return False, "No prenyl-derived side chain found"
    
    return True, "Contains quinone structure with prenyl-derived side chain"