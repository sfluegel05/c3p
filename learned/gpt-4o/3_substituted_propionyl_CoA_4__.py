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
    
    # More generalized CoA substructure pattern, targeting the pantetheine and core CoA motifs
    coa_core_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)C[C@H](O)C(C)(C)COP(O)(O)O")
    if not mol.HasSubstructMatch(coa_core_pattern):
        return False, "CoA backbone structure not found or incorrectly matched"
    
    # Broader thioester linkage pattern to accommodate variants
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)C")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage not adequately detected"
    
    # Identification of deprotonated phosphate groups with flexibility
    phosphate_pattern = Chem.MolFromSmarts("P([O-])(=O)([O-])O")
    matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(matches) < 2:
        return False, "Insufficient deprotonated phosphates indicating -4 charge"

    return True, "Molecule is classified as a 3-substituted propionyl-CoA(4-)"