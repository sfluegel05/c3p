"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: CHEBI:51625 nonclassic icosanoid

A nonclassic icosanoid is any biologically active signaling molecule made by oxygenation
of C20 fatty acids other than the classic icosanoids (the leukotrienes and the prostanoids).
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check carbon count (18-22)
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 18 or n_carbons > 22:
        return False, f"Carbon count ({n_carbons}) outside the expected range for nonclassic icosanoids"

    # Check for at least 2 oxygens
    n_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if n_oxygens < 2:
        return False, "Insufficient oxygenation for nonclassic icosanoid"

    # Check for specific oxygenated functional groups and unsaturation patterns
    epoxide_pattern = Chem.MolFromSmarts("[C@]1[C@@]2[C@@]1[C@]2")
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    conj_diene_pattern = Chem.MolFromSmarts("/C=C/C=C")

    epoxide_match = mol.HasSubstructMatch(epoxide_pattern)
    hydroxy_match = mol.HasSubstructMatch(hydroxy_pattern)
    conj_diene_match = mol.HasSubstructMatch(conj_diene_pattern)

    if not (epoxide_match or hydroxy_match or conj_diene_match):
        return False, "No characteristic oxygenation or unsaturation patterns found"

    # Check for specific substructures found in nonclassic icosanoids
    substructure_patterns = [
        Chem.MolFromSmarts("[C@@H]1[C@]2[C@@]1[C@]2"),  # epoxide adjacent to double bond
        Chem.MolFromSmarts("[C@H](O)[C@@H](O)"),  # vicinal diols
        Chem.MolFromSmarts("[C@H]([C@@H](O)O)O"),  # triol
        Chem.MolFromSmarts("C/C=C/C=C/C=C/C=C/C"),  # long conjugated system
        Chem.MolFromSmarts("C(=O)O"),  # carboxylic acid
    ]

    has_substructure = any(mol.HasSubstructMatch(pattern) for pattern in substructure_patterns)

    if not has_substructure:
        return False, "No characteristic substructures found for nonclassic icosanoid"

    return True, "Contains characteristic oxygenation, unsaturation, and substructures of nonclassic icosanoids"