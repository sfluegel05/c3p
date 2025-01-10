"""
Classifies: CHEBI:17984 acyl-CoA
"""
from rdkit import Chem

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is a thioester formed from the condensation of a carboxylic acid and coenzyme A.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Thioester linkage pattern: C(=O)S
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"
    
    # Adenine + Ribose + Diphosphate (part of the CoA backbone)
    ribose_adenine_pattern = Chem.MolFromSmarts("n1cncnc1C2OC(C(O)C2O)COP(=O)(O)OP(=O)(O)OC3COC3")
    if not mol.HasSubstructMatch(ribose_adenine_pattern):
        return False, "No adenine-ribose-diphosphate moiety found"
    
    # Pantetheine moiety including thiol
    pantetheine_pattern = Chem.MolFromSmarts("NCC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)O")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "No pantetheine moiety found"

    return True, "Contains thioester linkage and Coenzyme A structure"