"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Relaxed glycerol backbone with potential chirality at sn-2
    glycerol_pattern = Chem.MolFromSmarts("[C@H](O)C(O)CO") 
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with chirality at sn-2 found"

    # Phosphate group pattern accommodating more flexible phospho- motifs
    phosphate_pattern = Chem.MolFromSmarts("COP(O)(=O)OCCN")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No broad phosphoethanolamine group found"
    
    # Broadened ester linkage at sn-1
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage detected at sn-1"

    # Flexible acyl chain pattern
    acyl_chain_pattern = Chem.MolFromSmarts("C(=O)C[C,C]")
    if not mol.HasSubstructMatch(acyl_chain_pattern):
        return False, "No acyl chain detected"

    return True, "Molecule matches 1-acyl-sn-glycero-3-phosphoethanolamine structural criteria"