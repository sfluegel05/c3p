"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:15529 3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Check for CoA core structure
    adenine = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(adenine):
        return False, "No adenine moiety found"

    # Check for phosphate groups with specific patterns
    diphosphate = Chem.MolFromSmarts("[O-]P(=O)(O)OP(=O)([O-])O")
    terminal_phosphate = Chem.MolFromSmarts("OP([O-])([O-])=O")
    
    if not (mol.HasSubstructMatch(diphosphate) and mol.HasSubstructMatch(terminal_phosphate)):
        return False, "Missing required phosphate groups"

    # Count total negative charges - should be 4
    negative_charge_pattern = Chem.MolFromSmarts("[O-]")
    if len(mol.GetSubstructMatches(negative_charge_pattern)) != 4:
        return False, "Must have exactly 4 negative charges"

    # Check for thioester linkage
    thioester = Chem.MolFromSmarts("C(=O)SC")
    if not mol.HasSubstructMatch(thioester):
        return False, "No thioester linkage found"

    # Check for 3-hydroxy group on fatty acid chain
    # Account for both R and S stereochemistry
    hydroxy_pattern = Chem.MolFromSmarts("[C,c][CH1]([OH1])CC(=O)S")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3-hydroxy group found in correct position"

    # Check for pantetheine part
    pantetheine = Chem.MolFromSmarts("CC(C)(COP(=O))CO")
    if not mol.HasSubstructMatch(pantetheine):
        return False, "Missing pantetheine moiety"

    # Check for ribose
    ribose = Chem.MolFromSmarts("C1[CH1](O)[CH1](O)[CH1](O)C1")
    if not mol.HasSubstructMatch(ribose):
        return False, "Missing ribose moiety"

    # Verify fatty acid chain length (should be at least 4 carbons)
    fatty_acid = Chem.MolFromSmarts("CCCC")
    if not mol.HasSubstructMatch(fatty_acid):
        return False, "Fatty acid chain too short"

    return True, "Contains CoA(4-) moiety with 3-hydroxy fatty acid chain attached via thioester bond"