"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA results from the formal condensation of the thiol group 
    of coenzyme A with the carboxy group of any 3-hydroxy fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Coenzyme A pattern focused on adenine, ribose, phosphates, and pantetheine
    coa_pattern = Chem.MolFromSmarts("O[C@H]1CO[P@](O)(=O)O[C@@H]1COP(=O)(OC)COP(=O)(OC)COC(=O)NCCSC(=O)CC")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A structure not found"
    
    # 3-hydroxy fatty acyl pattern: 3-hydroxy group on a carbon chain, likely with a thioester linkage
    hydroxy_acyl_pattern = Chem.MolFromSmarts("C[C@H](O)CC(=O)S")
    if not mol.HasSubstructMatch(hydroxy_acyl_pattern):
        return False, "3-hydroxy group pattern not found"

    return True, "Contains structure consistent with a 3-hydroxy fatty acyl-CoA"