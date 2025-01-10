"""
Classifies: CHEBI:16337 phosphatidic acid
"""
from rdkit import Chem

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid is characterized by a glycerol backbone in which one hydroxy group
    is esterified with phosphoric acid and the other two are esterified with fatty acids.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycerol backbone pattern: O-CH2-CH(OH)-CH2-O
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Missing glycerol backbone pattern"

    # Check for ester linkages: C(=O)O
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups; need at least 2 for fatty acids"

    # Look for a phosphate group: P(=O)(O)(O)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Missing phosphate group"

    # Check total count of phosphorous, which should be 1 in PA
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 1:
        return False, "Phosphatidic acid should have exactly 1 phosphorus atom"

    return True, "Structure matches phosphatidic acid with glycerol backbone and appropriate esterifications"