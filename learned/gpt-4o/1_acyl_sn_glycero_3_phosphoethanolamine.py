"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone with stereo center
    glycerol_pattern = Chem.MolFromSmarts("[C@H](O)COC(=O)")  # Adjusted stereo recognition
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No correctly configured glycerol backbone found"
    
    # Check for presence of a phosphate group with potential ethanolamine
    phosphate_variants = ["COP(=O)(O)OCCN", "COP(=O)(O)OCC[NH3+]"]  # Added protonated variant
    phosphate_found = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in phosphate_variants)
    if not phosphate_found:
        return False, "No phosphate group with ethanolamine found"
    
    # Check for 1-O-acyl group (ester linkage)
    acyl_pattern = Chem.MolFromSmarts("COC(=O)")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No 1-O-acyl ester linkage found"
    
    return True, "Molecule fits all structural criteria for 1-acyl-sn-glycero-3-phosphoethanolamine"