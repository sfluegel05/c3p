"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
from rdkit import Chem

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol has a phosphatidyl group esterified to one of the hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Confirm inositol ring with 5 hydroxyl groups
    inositol_oh_pattern = Chem.MolFromSmarts("C1([C@@H]([C@H]([C@H]([C@@H]([C@H]1O)O)O)O)O)O") 
    if not mol.HasSubstructMatch(inositol_oh_pattern):
        return False, "No inositol ring with 5 hydroxyl groups found"

    # Phosphodiester linkage (glycerophosphate connected to inositol)
    phosphodiester_pattern = Chem.MolFromSmarts("O[P](=O)(O)[O][C@H]([O])[C@H]([O])[C@H]") 
    if not mol.HasSubstructMatch(phosphodiester_pattern):
        return False, "No glycerophosphate linkage to inositol found"

    # Recognition of esterified fatty acid chains connected to glycerol
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]C(=O)O")  # Example pattern for ester linkages in phosphatidyl group
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, f"Found {len(ester_matches)} ester groups, need at least 1 connected to phosphatidyl group"
    
    # Ensure connection of long aliphatic tail indicating fatty acid presence
    long_chain_pattern = Chem.MolFromSmarts("C~C~C~C~C~C")  # Generic long chain alkane pattern
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long aliphatic chain found indicating phosphatidyl groups"
    
    return True, "Consistent with a phosphatidylinositol: inositol esterified with glycerophosphate connected to phosphatidyl groups"