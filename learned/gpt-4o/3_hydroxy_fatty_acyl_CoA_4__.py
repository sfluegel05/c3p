"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhanced 3-hydroxy fatty acyl structure pattern
    hydroxy_fatty_acid_pattern = Chem.MolFromSmarts("[C,c][C,c][C@H](O)C(=O)C")
    if not mol.HasSubstructMatch(hydroxy_fatty_acid_pattern):
        return False, "No 3-hydroxy fatty acyl structure found"

    # Enhanced CoA pattern with comprehensive coverage
    coa_pattern = Chem.MolFromSmarts("NCC(=O)CCNCC(=O)SC[C,c]C")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA motif not present"

    # Correct phosphate patterns to ensure both are deprotonated
    deprotonated_phosphate_pattern = Chem.MolFromSmarts("P(=O)([O-])[O-]")
    phosphate_matches = mol.GetSubstructMatches(deprotonated_phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, "Not enough deprotonated phosphate groups found"

    return True, "Molecule classified as 3-hydroxy fatty acyl-CoA(4-)"