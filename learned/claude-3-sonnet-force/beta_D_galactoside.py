"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: CHEBI:28015 beta-D-galactoside
"""
from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    A beta-D-galactoside is any D-galactoside having beta-configuration at its anomeric centre.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-galactoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for beta-D-galactoside substructure
    beta_gal_patterns = [
        Chem.MolFromSmarts("[C@@H]1([C@H]([C@@H](O[C@]1(O)CO)O)O)O"),  # Explicit stereochemistry
        # Add more SMARTS patterns as needed
    ]
    
    for pattern in beta_gal_patterns:
        beta_gal_matches = mol.GetSubstructMatches(pattern)
        if beta_gal_matches:
            break
    else:
        return False, "No beta-D-galactoside substructure found"
    
    # Look for common substituents or modifications
    mod_patterns = [
        Chem.MolFromSmarts("[C@H](O)[C@H](CO)[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"),  # alpha-L-rhamnosyl
        Chem.MolFromSmarts("[C@H](O)[C@@H](CO)[C@H](O)[C@H](O)[C@H](O)"),  # alpha-L-arabinosyl
        Chem.MolFromSmarts("O=S(=O)(O)"),  # sulfate group
        Chem.MolFromSmarts("O=P(O)(O)O"),  # phosphate group
        # Add more patterns as needed
    ]
    
    for pattern in mod_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains beta-D-galactoside substructure with common modification"
    
    return True, "Contains unmodified beta-D-galactoside substructure"