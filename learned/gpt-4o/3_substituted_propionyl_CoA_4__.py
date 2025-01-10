"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem

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
    
    # Check for CoA substructure
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OCC1OC(C(O)C1O)N2C=NC3=C2N=CN=C3N")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA substructure not found"
    
    # Check for thioester linkage (general acyl group connected to CoA via a sulfur atom)
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage not found"
    
    # Check for any 3-substitution in propionyl group, which could be a carbon chain or other groups
    substitution_pattern = Chem.MolFromSmarts("CCC(=O)S")
    if not mol.HasSubstructMatch(substitution_pattern):
        return False, "No 3-substituted propionyl group found"
    
    # Check for the presence of deprotonated phosphate groups indicating a 4- charge
    phosphate_pattern = Chem.MolFromSmarts("P([O-])(=O)(O)OP([O-])(=O)O")
    matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(matches) < 2:
        return False, "Deprotonated phosphate groups indicating -4 charge not found"
    
    return True, "Molecule is a 3-substituted propionyl-CoA(4-)"