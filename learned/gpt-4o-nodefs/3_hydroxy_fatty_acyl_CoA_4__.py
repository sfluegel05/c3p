"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA_4_(smiles: str):
    """
    Classifies whether a chemical is a 3-hydroxy fatty acyl-CoA(4-).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Pattern to match hydroxyl group on fatty acyl chain at various positions with flexibility in stereochemistry and length
    hydroxyl_pattern = Chem.MolFromSmarts('[C@@H](O)CC(=O)SC*')  # Allow wildcard for variations in chain
    
    # Check for presence of the hydroxyl group at the 3rd position
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No 3-hydroxy fatty acyl chain found"
    
    # Comprehensive pattern for CoA moiety with flexibility for phosphate groups
    coa_pattern = Chem.MolFromSmarts('NC(=O)CNC(=O)C(COP([O-])(=O)O*)O[C@@H]1[C@H](O)C(COP([O-])([O-])=O)O[C@@H]2O[C@H](C(O)[C@@H](O)[C@H]2OP(*)*n2cnc3c(N)ncnc32)[C@H]1O')
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not fully matched"
    
    # Verify the carbon chain length to classify it as fatty acyl
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:  # Further relax the lower limit while ensuring it can capture longer hydrocarbon sequences
        return False, "Carbon chain too short for fatty acyl"
    
    return True, "Matches 3-hydroxy fatty acyl-CoA(4-) structure"