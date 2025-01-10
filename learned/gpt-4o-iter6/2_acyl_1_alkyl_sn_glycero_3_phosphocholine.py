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
    
    # Glycerol backbone with correct stereochemistry pattern
    sn_glycerol_pattern = Chem.MolFromSmarts("[C@H](COC)[O][CH]")
    if not mol.HasSubstructMatch(sn_glycerol_pattern):
        return False, "No sn-glycerol backbone with correct stereochemistry found"
    
    # Alkyl chain at sn-1 - ether linkage and long chain
    alkyl_chain_pattern = Chem.MolFromSmarts("COC[CX4]~[CX4]~[CX4]")
    if not mol.HasSubstructMatch(alkyl_chain_pattern):
        return False, "No alkyl chain at sn-1 position"
    
    # Acyl chain at sn-2 - ester linkage
    acyl_chain_pattern = Chem.MolFromSmarts("O[C@H](C)C(=O)~[CX4]")
    if not mol.HasSubstructMatch(acyl_chain_pattern):
        return False, "No acyl chain at sn-2 position"
    
    # Phosphocholine group: Phosphate and choline with quaternary nitrogen
    phosphocholine_pattern = Chem.MolFromSmarts("COP([O-])(=O)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"
    
    return True, "Contains sn-glycerol backbone with alkyl at sn-1, acyl at sn-2, and phosphocholine group"