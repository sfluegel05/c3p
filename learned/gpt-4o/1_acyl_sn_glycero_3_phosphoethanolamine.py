"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphoethanolamine group
    pe_pattern = Chem.MolFromSmarts("P(OCCN)(O)=O")
    if not mol.HasSubstructMatch(pe_pattern):
        return False, "No phosphoethanolamine group found"

    # Check for glycerol backbone with acyl linkage
    glycerol_pattern = Chem.MolFromSmarts("[C@H](O)COC(=O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with acyl linkage found"

    # Check for acyl group
    acyl_pattern = Chem.MolFromSmarts("C(=O)O[C@H]")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl group linked to chiral center found"

    return True, "Contains 1-acyl-sn-glycero-3-phosphoethanolamine structure"