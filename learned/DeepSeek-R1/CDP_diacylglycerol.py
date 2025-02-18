"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
"""
Classifies: CHEBI:17855 CDP-diacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    A CDP-diacylglycerol consists of a glycerol backbone with two acyl groups at positions 1 and 2,
    and a cytidine diphosphate (CDP) group at position 3.

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
    
    # Check for two ester groups (acyl chains)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"
    
    # Check glycerol backbone connected to diphosphate
    # Glycerol pattern: C-C-C with two esters and one diphosphate
    # SMARTS: [CH2]-[CH](-O-C(=O)-*)-[CH2]-O-P
    glycerol_diphosphate_pattern = Chem.MolFromSmarts("[CH2]-[CH](-O-C(=O)-*)-[CH2]-O-P")
    if not mol.HasSubstructMatch(glycerol_diphosphate_pattern):
        return False, "Glycerol backbone with diphosphate not found"
    
    # Check diphosphate bridge (P-O-P)
    diphosphate_pattern = Chem.MolFromSmarts("[PX4]-O-[PX4]")
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)
    if not diphosphate_matches:
        return False, "Diphosphate bridge not found"
    
    # Check cytidine moiety: cytosine + ribose connected to diphosphate
    # Cytosine pattern: pyrimidine ring with NH2 and C=O
    cytosine_pattern = Chem.MolFromSmarts("n1c(=O)[nH]c(N)cc1")
    if not mol.HasSubstructMatch(cytosine_pattern):
        return False, "Cytosine base not found"
    
    # Ribose pattern: five-membered ring with O and hydroxyl groups
    ribose_pattern = Chem.MolFromSmarts("[C@H]1O[C@H](CO)[C@H](O)[C@@H](O)[C@@H]1O")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "Ribose sugar not found"
    
    # Check connectivity between ribose and diphosphate
    # Combined pattern: ribose connected to diphosphate which is connected to glycerol
    # This SMARTS is complex and approximate; might need adjustment
    cdp_pattern = Chem.MolFromSmarts("[C@H]1O[C@H](COP(=O)(O)OP(=O)(O)O[C@H](COC(=O)*)C(OC(=O)*)CO)[C@H](O)[C@@H]1O")
    if not mol.HasSubstructMatch(cdp_pattern):
        return False, "CDP-glycerol structure not confirmed"
    
    # Optional: Check acyl chain length (e.g., at least 8 carbons each)
    # This is example-based and may vary
    # For each ester group, check the attached chain length
    # This part is commented out as chain length may vary
    # for ester in ester_matches:
    #     chain_length = 0
    #     # Traverse the acyl chain
    #     # Implementation depends on exact structure
    #     if chain_length < 8:
    #         return False, "Acyl chain too short"
    
    return True, "Contains glycerol with two acyl groups and CDP group"