"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
from rdkit import Chem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid must have two carboxylic acid groups and a hydroxy group
    on the alpha carbon to at least one of the carboxylic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-hydroxydicarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for carboxylic acid (-C(=O)O)
    carboxy_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    
    # Ensure we have exactly two carboxylic acid groups
    if len(carboxy_matches) < 2:
        return False, f"Expected at least 2 carboxylic acid groups but found {len(carboxy_matches)}"

    # Pattern for the alpha-hydroxy: HO-C-C(=O)O, accounting for chirality symbols
    alpha_hydroxy_pattern = Chem.MolFromSmarts("[C;H0,H1,H2](O)[C;H1,H2](C(=O)O)")
    alpha_oh_matches = mol.GetSubstructMatches(alpha_hydroxy_pattern)
    
    if not alpha_oh_matches:
        return False, "No hydroxy group found on the alpha carbon to a carboxylic acid group"

    return True, "Contains two carboxylic acid groups and a hydroxy group on the alpha carbon"