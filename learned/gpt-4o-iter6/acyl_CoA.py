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
    
    # Corrected pattern for thioester bond: C(=O)S
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"
    
    # Coenzyme A component: Pantetheine and adenine moieties
    pantetheine_pattern = Chem.MolFromSmarts("NCC(=O)CCNC(=O)[C@H](O)[C@H](C)[C@H](O)COP(O)O")
    adenine_moiety_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")  # Simple adenine structure
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "No pantetheine moiety found"
    if not mol.HasSubstructMatch(adenine_moiety_pattern):
        return False, "No adenine moiety found"

    return True, "Contains thioester linkage and Coenzyme A structure"