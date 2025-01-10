"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    Lysophosphatidic acids are characterized by a phosphate group, a glycerol backbone,
    and one fatty acid chain through ester linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a phosphate group
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for glycerol backbone - variations like non-chiral forms
    glycerol_pattern1 = Chem.MolFromSmarts("OCC(O)CO")
    glycerol_pattern2 = Chem.MolFromSmarts("[C@H](O)CO")
    if not (mol.HasSubstructMatch(glycerol_pattern1) or mol.HasSubstructMatch(glycerol_pattern2)):
        return False, "No appropriate glycerol backbone found"

    # Look for at least one ester linkage
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester linkage found for fatty acids"

    # Check for a sufficient number of carbons indicating a fatty acid
    c_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if c_count < 12:
        return False, "The molecule has too few carbons for a typical lysophosphatidic acid"

    return True, "Structure consistent with lysophosphatidic acid: contains phosphate group, glycerol backbone, and ester linkage for fatty acids"