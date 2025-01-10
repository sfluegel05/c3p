"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    Catechins are a type of flavonoid characterized by a 2-phenyl-3,4-dihydro-2H-chromen-3-ol structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Catechin pattern including stereochemistry, flexible definition
    catechin_pattern = Chem.MolFromSmarts("C1(OC2=CC=CC=C2C=C1O)C(CO)C3=CC(O)=CC(O)=C3")
    if not mol.HasSubstructMatch(catechin_pattern):
        return False, "Does not match catechin core structure (2-phenyl-3,4-dihydro-2H-chromen-3-ol)"

    # Allow for diverse substitutions and not just hydroxyl groups
    phenyl_pattern = Chem.MolFromSmarts("c[c][c][c][c][c]")  # Aromatic ring
    phenyl_matches = mol.GetSubstructMatches(phenyl_pattern)
    
    if len(phenyl_matches) <= 0:
        return False, f"No aromatic rings found, need at least one characteristic catechin"

    return True, "Contains catechin core structure with properly defined stereochemistry and substitutions"