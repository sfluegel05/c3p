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
    
    # Check for the core sn-glycerol backbone with invariant stereochemistry
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](O*)(COP([O-])(=O)OCC[N+](C)(C)C)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No SN-glycerol backbone with correct stereochemistry and phosphocholine group found"
   
    # Look for a long ether-linked alkyl chain from sn-1
    alkyl_ether_pattern = Chem.MolFromSmarts("COCCCCCCCCCCCCCCCC*")
    if not mol.HasSubstructMatch(alkyl_ether_pattern):
        return False, "No long alkyl ether chain detected at sn-1 position"
    
    # Check for ester-linked acyl chain at sn-2 with long chain attached
    acyl_ester_pattern = Chem.MolFromSmarts("OC([CX3](=O)CCCCCCCCCCCCC)*")
    if not mol.HasSubstructMatch(acyl_ester_pattern):
        return False, "No acyl chain at sn-2 position with ester linkage"
    
    # Confirm phosphocholine group with quaternary nitrogen
    phosphocholine_pattern = Chem.MolFromSmarts("P([O-])(=O)(OCC[N+](C)(C)C)")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group with quaternary nitrogen found"
    
    return True, "Contains sn-glycerol backbone with alkyl at sn-1, acyl at sn-2, and phosphocholine group"

# Test by passing examples of SMILES strings to the function