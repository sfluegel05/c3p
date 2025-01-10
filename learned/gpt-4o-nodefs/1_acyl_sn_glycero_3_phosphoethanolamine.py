"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    This is identified by a glycerol backbone, acyl group, and phosphoethanolamine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Glycerol backbone pattern (C-C-C with two oxygens)
    glycerol_pattern = Chem.MolFromSmarts("[C@H](O)[CH2]O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Acyl group - typically a long carbon chain connected to the glycerol backbone via ester bond
    acyl_pattern = Chem.MolFromSmarts("C(=O)O[C@H]")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl group found"
    
    # Phosphoethanolamine group - phosphate group with ethanolamine attached
    phospho_ethanolamine_pattern = Chem.MolFromSmarts("COP(O)(=O)OCCN")
    if not mol.HasSubstructMatch(phospho_ethanolamine_pattern):
        return False, "No phosphoethanolamine group found"

    # If all features are present, classify as 1-acyl-sn-glycero-3-phosphoethanolamine
    return True, "Contains glycerol backbone, acyl group, and phosphoethanolamine group"