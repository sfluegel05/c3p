"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    A lysophosphatidic acid has a glycerol backbone, one acyl chain, and a phosphate group,
    with no additional phospholipid head groups like choline.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone - allowing flexible stereochemistry
    # This encompasses several possible configurations to catch stereoisomers
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for the phosphate group; phosphoric acid monoester
    phosphate_patterns = [
        Chem.MolFromSmarts("OP(=O)(O)O"),  # Regular phosphate
    ]
    if not any(mol.HasSubstructMatch(p) for p in phosphate_patterns):
        return False, "No phosphate group found"

    # Look for monoacyl group (-O-C(=O)R) - checking for esters
    monoacyl_pattern = Chem.MolFromSmarts("C(=O)O")
    acyl_matches = mol.GetSubstructMatches(monoacyl_pattern)
    if len(acyl_matches) != 1:
        return False, f"Found {len(acyl_matches)} acyl groups, need exactly 1"

    # Ensure no head groups like choline are present
    choline_pattern = Chem.MolFromSmarts("OCC[N+](C)(C)C")
    if mol.HasSubstructMatch(choline_pattern):
        return False, "Choline head group found, not a lysophosphatidic acid"

    return True, "Contains glycerol backbone, one acyl group, and a phosphate group without additional head groups"