"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    A catechin generally has a 2-phenyl-3,4-dihydro-2H-chromen-3-ol structure which may include various substitutions.

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
    
    # Define generic catechin pattern allowing for core flexibility and substitutions
    # The general structure is C1(OC2=C(O)C=CC=C2C=C1)C(O)C3=CC(O)=CC(O)=C3
    catechin_patterns = [
        Chem.MolFromSmarts("C1(OC2=CC=CC=C2C=C1)C(CO)C3=CC(O)=CC(O)=C3"),      # Base pattern
        Chem.MolFromSmarts("C1(OC2=C(C=CC=C2)C=C1)C(O)C3=CC(O)=CC(O)=C3"),    # Flexibly accommodating rearrangements
        Chem.MolFromSmarts("C1(OC2=CC=CC=C2C=C1)C(CO)C3=C(O)C(O)=C(O)C=C3"),  # Including additional hydroxyls
    ]

    if not any(mol.HasSubstructMatch(pattern) for pattern in catechin_patterns):
        return False, "Does not match catechin core structure (2-phenyl-3,4-dihydro-2H-chromen-3-ol)"
    
    return True, "Contains catechin core structure with permissible stereochemical and substitution variations"