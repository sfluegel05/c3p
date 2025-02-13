"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
from rdkit import Chem

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.
    A phosphatidylethanolamine is a glycerophospholipid where two fatty acids are esterified to a glycerol backbone,
    and a phosphatidyl group is esterified to ethanolamine.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylethanolamine, False otherwise.
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone including the phosphate group
    glycerol_backbone_pattern = Chem.MolFromSmarts("[C@H](O)(CO[P](=O)(O)[O-])")
    if not mol.HasSubstructMatch(glycerol_backbone_pattern):
        return False, "No glycerol-like backbone with phosphate group attached found"

    # Check for two esterified fatty acids
    ester_pattern = Chem.MolFromSmarts("COC(=O)C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester linkages, need at least 2"

    # Check for phosphoethanolamine group
    phosphoethanolamine_pattern = Chem.MolFromSmarts("OP(O)(OCCN)O")
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "No phosphoethanolamine group found"
    
    # Successfully matched all criteria
    return True, "Contains glycerol backbone with 2 esterified fatty acids and a phosphoethanolamine group"

# Example usage:
# is_phosphatidylethanolamine("CCCCCCCCCCCCCCCCCC(=O)O[C@H](COC(=O)CCCCCCC\C=C/CCCCCCCC)COP(O)(=O)OCCN")