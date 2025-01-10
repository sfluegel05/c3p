"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    Key features:
    - Contains Coenzyme A moiety
    - Contains a thioester linkage
    - Contains a 3-oxo group in the fatty acyl chain

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Update the CoA moiety SMARTS pattern to incorporate more structural elements
    # Consider all features—ribose, phosphate linkages, adenine, pantetheine moiety
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@@H](O)C(C)(C)COP(O)(=O)O[C@@H]1CO[C@H](O)[C@H]1O[P,S](=O)(O)OP(=O)(O)C[C@H]2O[C@H]([C@H](O)[C@H]2OP(=O)(O)O)n3cnc4c3ncnc4N")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"
    
    # Thioester linkage
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)C")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"
    
    # 3-oxo (ketone) group pattern at the third carbon of acyl chain
    oxo_pattern = Chem.MolFromSmarts("CCC(=O)C")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 3-oxo group found in the fatty acyl chain"
    
    return True, "Contains CoA moiety, thioester linkage, and 3-oxo group in the fatty acyl chain"