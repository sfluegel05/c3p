"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: catechols (compounds containing an o-diphenol component)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catechols(smiles: str):
    """
    Determines if a molecule contains a catechol group (o-diphenol).
    A catechol has two hydroxyl groups (-OH) in ortho position on a benzene ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a catechol group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for catechol group:
    # Two adjacent carbons on aromatic ring, each with hydroxyl group
    catechol_pattern = Chem.MolFromSmarts("c1c(O)c(O)cccc1")
    
    # Alternative pattern that also matches other arrangements
    alt_pattern = Chem.MolFromSmarts("c([OH1])c([OH1])")
    
    # Check for matches
    if mol.HasSubstructMatch(catechol_pattern) or mol.HasSubstructMatch(alt_pattern):
        # Get number of matches to include in reason
        matches = len(mol.GetSubstructMatches(catechol_pattern))
        alt_matches = len(mol.GetSubstructMatches(alt_pattern))
        total_matches = max(matches, alt_matches)
        
        if total_matches == 1:
            return True, "Contains one catechol group (o-diphenol)"
        else:
            return True, f"Contains {total_matches} catechol groups"
    
    # Check if molecule has any aromatic rings
    if not any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "No aromatic rings found"
        
    # Check if molecule has any hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH1]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl groups found"
        
    return False, "No catechol pattern (adjacent hydroxyl groups on benzene ring) found"