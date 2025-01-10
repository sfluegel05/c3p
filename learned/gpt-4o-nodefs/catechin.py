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
    # The general catechin structure is given as: 2-phenyl-3,4-dihydro-2H-chromen-3-ol
    catechin_patterns = [
        Chem.MolFromSmarts("C1(C=C2c3cc(O)ccc3OC2CC1)"),             # Basic catechin structure
        Chem.MolFromSmarts("C1(O[C@H](cc2ccccc2)[C@@H](O)CO1)")       # Configuration with hydroxyls and flexible stereochemistry
    ]

    # Check each pattern for substructure match
    if not any(mol.HasSubstructMatch(pattern) for pattern in catechin_patterns):
        return False, "Does not match catechin core structure (2-phenyl-3,4-dihydro-2H-chromen-3-ol)"
    
    return True, "Contains catechin core structure with permissible stereochemical and substitution variations"