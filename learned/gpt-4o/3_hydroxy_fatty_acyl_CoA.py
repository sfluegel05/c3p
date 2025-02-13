"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA is the result of condensation of the thiol group of coenzyme A
    with the carboxy group of a 3-hydroxy fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Coenzyme A pattern including ribose, adenine, and phosphates
    coa_pattern = Chem.MolFromSmarts("NC(=O)C(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)C(O)=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A structure not found"
    
    # Pattern for 3-hydroxy fatty acyl: 3-hydroxy group bonded to a carbon chain, linked by a thioester bond
    hydroxy_acyl_pattern = Chem.MolFromSmarts("C[C@@H](O)C(=O)SCC")
    if not mol.HasSubstructMatch(hydroxy_acyl_pattern):
        return False, "3-hydroxy group with correct linkage not found"

    return True, "Contains structure consistent with a 3-hydroxy fatty acyl-CoA"