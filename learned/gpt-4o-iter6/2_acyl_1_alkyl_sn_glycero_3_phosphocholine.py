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
    
    # Correct stereochemistry of sn-glycerol backbone
    # This pattern assumes variable ether linkages with stereocenter
    glycerol_pattern = Chem.MolFromSmarts("O[C@H](COP([O-])(=O)OCC[N+](C)(C)C)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No SN-glycerol backbone with correct stereochemistry and phosphocholine group found"
   
    # Alkyl chain at sn-1 - commonly via ether bond, check for ether link but allow variations
    alkyl_chain_pattern = Chem.MolFromSmarts("CO[CH2][CX4,H]")
    if not mol.HasSubstructMatch(alkyl_chain_pattern):
        return False, "No alkyl ether chain at sn-1 position detected"
    
    # Acyl chain at sn-2 - typically an ester linkage
    acyl_chain_pattern = Chem.MolFromSmarts("C(O[CX4])(=O)~[CX4,CH]")
    if not mol.HasSubstructMatch(acyl_chain_pattern):
        return False, "No acyl chain at sn-2 position with ester linkage"
    
    # Check for existence of quaternary ammonium ion pattern in phosphocholine
    phosphocholine_pattern = Chem.MolFromSmarts("P([O-])(=O)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group with quaternary nitrogen found"
    
    return True, "Contains sn-glycerol backbone with alkyl at sn-1, acyl at sn-2, and phosphocholine group"