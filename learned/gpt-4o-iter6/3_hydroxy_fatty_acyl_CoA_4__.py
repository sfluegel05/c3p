"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA_4_(smiles: str):
    """
    Classifies a molecule as a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.
    
    The criteria focus on the presence of a 3-hydroxy group, a Coenzyme A moiety, 
    and deprotonated phosphate groups. It also considers the long fatty acid chain.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a 3-hydroxy fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Check for a 3-hydroxy group with correct positioning
    three_hydroxy_pattern = Chem.MolFromSmarts("[CX4,c][C@H](O)[CX3](=O)")
    if not mol.HasSubstructMatch(three_hydroxy_pattern):
        return False, "No correctly positioned 3-hydroxy group found on fatty acid chain"

    # Coenzyme A moiety pattern
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)C(C)(C)COP(=O)([O-])O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A moiety with thioester linkage not found"

    # Look for deprotonated phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([O-])[O-]")
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_matches < 2:
        return False, f"Found {phosphate_matches} deprotonated phosphate groups, need at least 2"

    # Fatty acid chain pattern (accounting for variability in chain length and saturation)
    chain_length_pattern = Chem.MolFromSmarts("[C@]([O])(C=O)SCCN")  # Better representation of 3-hydroxy-palmitoyl structural diversification
    if not mol.HasSubstructMatch(chain_length_pattern):
        return False, "Insufficient long and structurally correct fatty acyl chain detected"

    return True, "The molecule is a 3-hydroxy fatty acyl-CoA(4-)"