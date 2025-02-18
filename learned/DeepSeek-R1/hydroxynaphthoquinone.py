"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: CHEBI:82525 hydroxynaphthoquinone
"""
from rdkit import Chem
from rdkit.Chem import Mol

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    A hydroxynaphthoquinone is a naphthoquinone with at least one hydroxyl group substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxynaphthoquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for naphthoquinone structure: two ketone groups on naphthalene
    quinone_pattern = Chem.MolFromSmarts("c1c2ccccc2c(c1=O)=O")  # Naphthoquinone pattern
    if not mol.HasSubstructMatch(quinone_pattern):
        return False, "No naphthoquinone moiety found"

    # Check for at least one hydroxyl group (-OH) attached to the aromatic rings
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]-c")  # Hydroxyl attached to aromatic carbon
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing hydroxyl substituent on aromatic ring"

    return True, "Contains naphthoquinone moiety with at least one hydroxyl group"