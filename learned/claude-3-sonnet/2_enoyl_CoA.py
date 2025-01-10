"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA is an unsaturated fatty acyl-CoA in which the S-acyl group 
    contains a double bond between positions 2 and 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA backbone pattern
    # Look for adenine + ribose phosphate + pantetheine parts
    coa_pattern = Chem.MolFromSmarts("[nX2r6:1]1[cX3r6:2][nX2r6:3][cX3r6:4][cX3r6:5][nX2r6:6]1")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found (missing adenine)"
    
    # Look for thioester group (-C(=O)S-)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatches(thioester_pattern):
        return False, "No thioester group found"

    # Look for double bond in position 2-3 relative to the thioester
    # Pattern: -C(=O)S-CH2-CH2-NH-C(=O)-CH2-CH2-NH-C(=O)- (pantetheine part)
    enoyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2][CH2][CH2]")
    if not mol.HasSubstructMatch(enoyl_pattern):
        return False, "Missing pantetheine part"

    # Look for the characteristic C=C-C(=O)S pattern of 2-enoyl-CoA
    alpha_beta_pattern = Chem.MolFromSmarts("[CX4,H;!$(C=O)]-[CX3]=[CX3]-C(=O)[SX2]")
    if not mol.HasSubstructMatch(alpha_beta_pattern):
        return False, "No double bond between positions 2 and 3 relative to thioester"

    # Additional check for phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=[OX1])[OX2H,OX1-]")
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_matches < 3:
        return False, "Missing phosphate groups characteristic of CoA"

    return True, "Contains CoA moiety and unsaturated fatty acyl group with double bond between positions 2 and 3"