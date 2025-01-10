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
    
    # Pattern for identifying 3-hydroxy group in fatty acyl chain
    hydroxyl_pattern = Chem.MolFromSmarts('[CX4,CX3][C@@H](O)C(=O)SC')
    
    # Check for presence of the hydroxyl group at the 3rd position
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No 3-hydroxy fatty acyl chain found"

    # Pattern for Coenzyme A structure with the nucleotide section and phosphate groups
    coa_pattern = Chem.MolFromSmarts('NC(=O)CCNC(=O)O[C@@H]1[C@H](O)[C@@H]([C@@H]1OP(=O)([O-])O)O[C@@H]1[C@H](O)[C@H](COP(=O)([O-])O)O1')
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not fully matched"
    
    # Verify the carbon chain length for the fatty acyl component
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Carbon chain too short for fatty acyl"

    return True, "Matches 3-hydroxy fatty acyl-CoA(4-) structure"