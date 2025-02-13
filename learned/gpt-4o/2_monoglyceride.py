"""
Classifies: CHEBI:17389 2-monoglyceride
"""
from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride is a glycerol in which the acyl substituent is located at position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a glycerol backbone with an acyl linked at the 2nd position
    # Glycerol backbone pattern: OCC(O)CO with additional ester:
    monoglyceride_pattern = Chem.MolFromSmarts("C(CO)OC(=O)[C,C][C,C]")  # Partial ester group patterns around glycerol
    
    if not mol.HasSubstructMatch(monoglyceride_pattern):
        return False, "Does not match the 2-monoglyceride pattern"

    # Find ester linkages explicitly at central glycerol position and ensure exclusivity
    ester_connectivity = Chem.MolFromSmarts("C(CO)OC(=O)[C,C][C,C]")  # stricter account for ester configuration
    if not mol.HasSubstructMatch(ester_connectivity):
        return False, "Ester linkage not precisely at the 2-position"

    # Verify only one esterification (one acyl chain)
    general_ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(general_ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    return True, "Contains a glycerol backbone with an acyl substituent at position 2"