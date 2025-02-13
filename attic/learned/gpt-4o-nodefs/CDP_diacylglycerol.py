"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    A CDP-diacylglycerol has a glycerol backbone, two ester-linked fatty acid chains and a cytidine diphosphate group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a CDP-diacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for glycerol backbone, ester groups, and cytidine diphosphate
    glycerol_backbone_pattern = Chem.MolFromSmarts("OCC(O)CO")
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    cytidine_pattern = Chem.MolFromSmarts("N2C=CC(=NC2=O)N")
    diphosphate_pattern = Chem.MolFromSmarts("P(OP(=O)(O)O)(=O)")
    
    # Check for glycerol backbone
    if not mol.HasSubstructMatch(glycerol_backbone_pattern):
        return False, "No glycerol backbone found"
        
    # Check for at least two ester groups (two fatty acid chains)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"
    
    # Check for cytidine structure
    if not mol.HasSubstructMatch(cytidine_pattern):
        return False, "No cytidine moiety found"
        
    # Check for diphosphate group
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "No cytidine diphosphate group found"
    
    return True, "Contains glycerol backbone with two fatty acid chains and cytidine diphosphate group"