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

    # Check for adenine core
    adenine = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(adenine):
        return False, "No adenine moiety found"

    # Check for phosphate groups - more flexible pattern
    phosphate_pattern = Chem.MolFromSmarts("OP([O-])(=O)OP([O-])(=O)O")
    terminal_phosphate = Chem.MolFromSmarts("OP([O-])([O-])=O")
    
    if not (mol.HasSubstructMatch(phosphate_pattern) and mol.HasSubstructMatch(terminal_phosphate)):
        return False, "Missing required phosphate groups"

    # Count negative charges - should be 4
    negative_charges = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[O-]")))
    if negative_charges != 4:
        return False, f"Found {negative_charges} negative charges, need exactly 4"

    # Check for thioester linkage
    thioester = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(thioester):
        return False, "Missing or incorrect thioester-pantetheine linkage"

    # Check for 3-hydroxy group on fatty acid chain
    # More flexible pattern accounting for both R and S stereochemistry
    hydroxy_pattern = Chem.MolFromSmarts("[C,c][CH1]([OH1])CC(=O)S")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3-hydroxy group found in correct position"

    # Check for ribose with correct connectivity
    ribose = Chem.MolFromSmarts("[C@H]1O[C@H](COP)[C@H](O)[C@@H](OP)[C@H]1O")
    if not mol.HasSubstructMatch(ribose):
        return False, "Missing or incorrect ribose moiety"

    # Check for pantetheine part with more complete pattern
    pantetheine = Chem.MolFromSmarts("NC(=O)CCNC(=O)[CH](O)C(C)(C)COP")
    if not mol.HasSubstructMatch(pantetheine):
        return False, "Missing or incorrect pantetheine moiety"

    # Verify fatty acid chain (allowing for saturated and unsaturated chains)
    fatty_chain = Chem.MolFromSmarts("C[C,$(C=C)]C[C,$(C=C)]C")
    if not mol.HasSubstructMatch(fatty_chain):
        return False, "No suitable fatty acid chain found"

    # Additional check for complete CoA structure
    coa_backbone = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)C(O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP([O-])([O-])=O")
    if not mol.HasSubstructMatch(coa_backbone):
        return False, "Incomplete or incorrect CoA backbone structure"

    return True, "Contains complete 3-hydroxy fatty acyl-CoA(4-) structure with correct phosphate charges"