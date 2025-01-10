"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem

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
    
    # Recognize the glycerol backbone with chiral center at sn-2
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    chiral_center_at_sn2 = any(center[1] == 'R' for center in chiral_centers if center[0] == 1)  # Ideally index should match sn-2 position
    if not chiral_center_at_sn2:
        return False, "No glycerol backbone with (R)-configuration at sn-2 position found"

    # Check for the presence of glycerol-like structure
    glycerol_pattern = Chem.MolFromSmarts("[C@H](O)[C@H](O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol-like backbone found"

    # Confirm the attached phosphoethanolamine group
    phosphate_ethanolamine_pattern = Chem.MolFromSmarts("P(=O)(O)OCCN")
    if not mol.HasSubstructMatch(phosphate_ethanolamine_pattern):
        return False, "No phosphoethanolamine group found"
    
    # Detect ester bond location strictly at sn-1 position
    ester_sn1_pattern = Chem.MolFromSmarts("O[C@H](CO)C(=O)")  
    if not mol.HasSubstructMatch(ester_sn1_pattern):
        return False, "No ester linkage at sn-1 position detected"
    
    return True, "Molecule matches 1-acyl-sn-glycero-3-phosphoethanolamine structural criteria"