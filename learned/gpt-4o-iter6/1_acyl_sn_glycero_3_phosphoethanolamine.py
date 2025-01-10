"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    This class is characterized as a 1-O-acylglycerophosphoethanolamine with (R)-configuration at sn-2.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Recognize the glycerol backbone with chiral center at sn-2 (R-configuration)
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](O)CO[C@H](COP(=O)(O)OCCN)O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with right chirality at sn-2 found"

    # Search for a phosphate group connected to ethanolamine
    phosphoethanolamine_pattern = Chem.MolFromSmarts("COP(=O)(O)OCCN")
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "No phosphoethanolamine group found"
    
    # Look for ester linkage at sn-1
    ester_linkage_pattern = Chem.MolFromSmarts("C(=O)O[C@H]")
    if not mol.HasSubstructMatch(ester_linkage_pattern):
        return False, "No ester linkage detected at sn-1"
    
    # Identify any fatty acyl chain
    acyl_chain_pattern = Chem.MolFromSmarts("C(=O)C")
    if not mol.HasSubstructMatch(acyl_chain_pattern):
        return False, "No acyl chain detected"
    
    return True, "Molecule matches 1-acyl-sn-glycero-3-phosphoethanolamine structural criteria"