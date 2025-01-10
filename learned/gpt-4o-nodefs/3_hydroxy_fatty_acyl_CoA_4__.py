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
        
    # Pattern to match 3-hydroxy group on fatty acyl chain at various positions (R is wildcard for hydrocarbon chain)
    hydroxyl_pattern = Chem.MolFromSmarts('[C@@H](O)CC(=O)SC')  # Pattern includes stereochemistry and thioester link
    
    # Check for presence of the hydroxyl group at the 3rd position
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No 3-hydroxy fatty acyl chain found"
    
    # Comprehensive pattern for CoA moiety to include adenosine, phosphates, and thiol
    coa_pattern = Chem.MolFromSmarts('NC(=O)CNC(=O)C(COP(O)([O-])=O)O[C@@H]1[C@H](O)C(COP(O)(=O)O[C@H]2O[C@H](C(O)[C@@H](O)[C@H]2OP([O-])([O-])=O)n2cnc3c(N)ncnc32)[C@H]1O')
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not fully matched"
    
    # Check for sufficient length of carbon chain to be considered fatty acyl
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 14:  # Increasing possible lower limit to capture more examples
        return False, "Carbon chain too short for fatty acyl"
    
    return True, "Matches 3-hydroxy fatty acyl-CoA(4-) structure"