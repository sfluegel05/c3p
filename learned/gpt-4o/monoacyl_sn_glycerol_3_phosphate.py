"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    This class can have a single acyl group at either position 1 or position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacyl-sn-glycerol 3-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone with a phosphodiester
    glycerol_phosphate_pattern = Chem.MolFromSmarts("[O;H1:3][C@@H]([O:1])COP(O)(O)=O")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol backbone with phosphate group found"
    
    # Look for the presence of one ester linkage (-C(=O)OC-) on the glycerol backbone
    ester_pattern = Chem.MolFromSmarts("C(=O)OC[C@@H](COP(O)(O)=O)O")
    ester_match_1 = mol.HasSubstructMatch(ester_pattern)
    
    # Check alternative acyl attachment on the glycerol - generalized as position-1 or 2 attachment
    alt_ester_pattern = Chem.MolFromSmarts("C(=O)O[C@@H](CO)COP(O)(O)=O")
    ester_match_2 = mol.HasSubstructMatch(alt_ester_pattern)
    
    # Confirm the presence of only one ester linkage corresponding to only one acyl group
    if (ester_match_1 and ester_match_2) or (not ester_match_1 and not ester_match_2):
        return False, "Incorrect number of ester linkages: Expected exact one"
    
    return True, "Contains glycerol backbone with one acyl group and a phosphate group"