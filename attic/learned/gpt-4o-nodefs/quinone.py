"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    A quinone typically contains a fully conjugated cyclic dione structure,
    usually with carbonyl groups in 1,4 or 1,2 positions on an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for quinone pattern, e.g., 1,4-benzoquinone
    quinone_pattern = Chem.MolFromSmarts("O=C1C=CC(=O)C=C1")  # 1,4-benzoquinone pattern
    naphthoquinone_pattern = Chem.MolFromSmarts("O=C1C=CC2=CC=CC(=O)C2=C1")  # Naphthoquinone pattern
    
    if mol.HasSubstructMatch(quinone_pattern):
        return True, "Contains 1,4-benzoquinone substructure"
    elif mol.HasSubstructMatch(naphthoquinone_pattern):
        return True, "Contains naphthoquinone substructure"

    # Check generally for two carbonyl groups on an aromatic ring
    ketone_pattern = Chem.MolFromSmarts("C=O")
    aromatic_ring = Chem.MolFromSmarts("a1aaaaa1")  # basic aromatic ring pattern

    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    aromatic_matches = mol.GetSubstructMatches(aromatic_ring)
    
    if len(ketone_matches) >= 2 and any(len(aromatic) >= 6 for aromatic in aromatic_matches):
        return True, "Contains at least two ketone groups on an aromatic ring"

    return False, "No typical quinone structure matching pattern found"