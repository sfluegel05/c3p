"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is defined as a six-carbon monosaccharide with an aldehyde at position 1 (aldohexose)
    or a ketone group at position 2 (ketohexose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.HybridizationType.SP3)
    
    # Check if there are exactly six sp3-hybridized carbon atoms
    if c_count != 6:
        return False, f"Contains {c_count} sp3-hybridized carbon atoms, but a hexose requires exactly 6"

    # Define more specific patterns for aldehyde or ketone groups
    aldehyde_pattern = Chem.MolFromSmarts("[CH]=O")  # Aldehyde group
    ketone_pattern = Chem.MolFromSmarts("[C](=O)[C]")  # Ketone group in context

    # Check for an aldehyde at position 1 (linear form)
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Structure matches an aldohexose (aldehyde group present)"
    
    # Check for a ketone at position 2
    if mol.HasSubstructMatch(ketone_pattern):
        return True, "Structure matches a ketohexose (ketone group present)"

    # Enhanced cyclic form patterns (pyranose or furanose)
    furanose_pattern = Chem.MolFromSmarts("C1OC[C@H](O)[C@@H]1O")
    pyranose_pattern = Chem.MolFromSmarts("C1OC[C@H](O)[C@@H](O)C1")

    if mol.HasSubstructMatch(furanose_pattern) or mol.HasSubstructMatch(pyranose_pattern):
        return True, "Structure matches a cyclic furanose or pyranose form"

    return False, "Does not contain hexose-defining functional groups (aldehyde or ketone)"