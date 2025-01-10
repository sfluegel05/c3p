"""
Classifies: CHEBI:33447 phospho sugar
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar contains a sugar moiety (often as a furanose or pyranose)
    with one or more attached phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible phosphate group pattern with substructure query
    phosphate_patterns = [
        Chem.MolFromSmarts("[OX2]P(=O)([OX1])[OX1]"),  # General phosphate
        Chem.MolFromSmarts("[OX2]P(=O)([OX1])[OX1-]"),  # Deprotonated version
    ]

    # Check for at least one phosphate group pattern match
    has_phosphate = any(mol.HasSubstructMatch(pat) for pat in phosphate_patterns)
    if not has_phosphate:
        return False, "No phosphate group found"

    # Define sugar ring patterns (furanose or pyranose or any carbohydrate ring)
    sugar_ring_patterns = [
        Chem.MolFromSmarts("C1OC(CO)C(O)C1"),  # Pyranose
        Chem.MolFromSmarts("C1OC(C)C1"),       # Furanose
        Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1"),  # Open sugar structures
        Chem.MolFromSmarts("O=C[C@H](O)[C@@H](O)C=O"),  # Linear carbohydrate fragments
    ]

    # Checking for the presence of any sugar-related structure
    has_sugar_structure = any(mol.HasSubstructMatch(pat) for pat in sugar_ring_patterns)
    if not has_sugar_structure:
        return False, "No sugar ring structure detected"

    # If both sugar and phosphate are present, classify as phospho sugar
    return True, "Contains a phosphate group attached to a sugar moiety"