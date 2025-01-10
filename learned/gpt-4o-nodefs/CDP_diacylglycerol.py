"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
from rdkit import Chem

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
    
    # Refine SMARTS pattern for glycerol backbone with potential stereochemistry
    glycerol_backbone_pattern = Chem.MolFromSmarts("[C@H](O)[C@H](O)CO")
    
    # Adjust cytidine detection to account for potential flexibility in structure
    cytidine_pattern = Chem.MolFromSmarts("n1c(C)c[nH]c1=O")
    
    # Attempt a refined global pattern for phosphatidyl components including diphosphate
    diphosphate_pattern = Chem.MolFromSmarts("P(O[PH](O[CX4])[CX4])O")
    
    # Check for glycerol backbone
    if not mol.HasSubstructMatch(glycerol_backbone_pattern):
        return False, "No glycerol backbone found"
        
    # Check for ester linkages, looking for two
    ester_pattern = Chem.MolFromSmarts("C(=O)O[CH2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"
    
    # Check for cytidine moiety structure
    if not mol.HasSubstructMatch(cytidine_pattern):
        return False, "No cytidine moiety found"
        
    # Check for presence of a diphosphate group
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "No cytidine diphosphate group found"
    
    return True, "Contains glycerol backbone with two ester-linked fatty acid chains and cytidine diphosphate group"