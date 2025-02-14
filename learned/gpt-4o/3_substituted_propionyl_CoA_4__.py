"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-substituted propionyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Adjusted CoA substructure pattern with flexibility in the nucleotide part
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OCC1OC(C(O)C1O)N2C=NC3=C2N=CN=C3N")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA substructure not found or incorrectly matched due to specificity"
    
    # Thioester linkage pattern that allows for variations in the acyl chain
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage not found, or pattern is too specific"
    
    # Check for potential 3-substituted pattern, focusing on general acyl properties
    substitution_pattern = Chem.MolFromSmarts("CC(=O)S")
    if not mol.HasSubstructMatch(substitution_pattern):
        return False, "No appropriate 3-substituted propionyl-like group identified"
    
    # Check for deprotonated phosphate groups, indicating a 4- charge state
    phosphate_pattern = Chem.MolFromSmarts("P([O-])(=O)(O)OP([O-])(=O)O")
    matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(matches) < 2:
        return False, "Insufficient deprotonated phosphate groups to indicate -4 charge"

    return True, "Molecule is classified as a 3-substituted propionyl-CoA(4-)"