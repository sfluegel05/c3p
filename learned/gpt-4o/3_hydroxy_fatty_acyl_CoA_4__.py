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

    # 3-hydroxy fatty acid structure pattern, adjusted for potential stereoisomerism
    hydroxy_fatty_acid_pattern = Chem.MolFromSmarts("[C,c][C,c][C@@H,O](O)C(=O)[C,c][C,c]")
    if not mol.HasSubstructMatch(hydroxy_fatty_acid_pattern):
        return False, "No 3-hydroxy fatty acyl structure found"

    # Comprehensive CoA pattern covering all relevant functional groups
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNCC(=O)SC[C,c][C@@H](O)CoAP")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA motif not present"

    # Verify deprotonation with a more comprehensive pattern
    deprotonated_phosphate_pattern = Chem.MolFromSmarts("P(=O)([O-])[O-]")
    if len(mol.GetSubstructMatches(deprotonated_phosphate_pattern)) < 2:
        return False, "Not enough deprotonated phosphate groups found"

    return True, "Molecule classified as 3-hydroxy fatty acyl-CoA(4-)"