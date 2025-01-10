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
    pe_pattern = Chem.MolFromSmarts("COP(O)(=O)OCCN")
    if not mol.HasSubstructMatch(pe_pattern):
        return False, "No phosphoethanolamine group found"

    # Check for glycerol backbone with an acyl linkage
    glycerol_acyl_pattern = Chem.MolFromSmarts("[C@@H](O)COP(O)(=O)OCCN")  # Stereochemistry specified
    if not mol.HasSubstructMatch(glycerol_acyl_pattern):
        return False, "No glycerol backbone with chiral center found"
    
    # Check specifically for acyl linkage indicating ester
    acyl_linkage_pattern = Chem.MolFromSmarts("O[C@@H](CO)C(=O)C")
    if not mol.HasSubstructMatch(acyl_linkage_pattern):
        return False, "No acyl group properly linked to glycerol backbone found"

    return True, "Contains 1-acyl-sn-glycero-3-phosphoethanolamine structure"