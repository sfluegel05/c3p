"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for glycerol backbone with sn stereochemistry pattern
    sn_glycerol_pattern = Chem.MolFromSmarts("[C@H](CO*)(O*)")
    if not mol.HasSubstructMatch(sn_glycerol_pattern):
        return False, "No sn-glycerol backbone found"
    
    # Check for Alkyl chain at sn-1 position
    alkyl_chain_pattern = Chem.MolFromSmarts("COC")
    if not mol.HasSubstructMatch(alkyl_chain_pattern):
        return False, "No alkyl chain at sn-1 position"
    
    # Check for Acyl chain at sn-2 position
    acyl_chain_pattern = Chem.MolFromSmarts("[C](O*)(=O)*")
    if not mol.HasSubstructMatch(acyl_chain_pattern):
        return False, "No acyl chain at sn-2 position"
    
    # Check for phosphate and choline group
    phosphocholine_pattern = Chem.MolFromSmarts("COP([O-])(=O)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphate and choline group"
    
    return True, "Contains sn-glycerol backbone with alkyl at sn-1, acyl at sn-2, and phosphocholine group"