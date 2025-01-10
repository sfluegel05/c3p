"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    A hydroxynaphthoquinone is a naphthoquinone moiety with at least one hydroxy group substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxynaphthoquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expanded SMARTS pattern for naphthoquinone with flexibility
    naphthoquinone_patterns = [
        Chem.MolFromSmarts('O=C1C=CC(=O)c2ccccc12'),  # core naphthoquinone
        Chem.MolFromSmarts('O=C1CC(=O)c2cccc(O)c12'), # alternative patterns
        Chem.MolFromSmarts('c1cc(C(=O)c2ccccc2C1=O)O'), # presence with hydroxyl
    ]

    # Check if any pattern matches
    if not any(mol.HasSubstructMatch(pattern) for pattern in naphthoquinone_patterns if pattern):
        return False, "No valid naphthoquinone core found"

    # Verify presence of hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts('[OH]')
    if mol.HasSubstructMatch(hydroxy_pattern):
        return True, "Contains naphthoquinone core with hydroxy group substitution"

    return False, "Naphthoquinone core present but no hydroxy substitution found"