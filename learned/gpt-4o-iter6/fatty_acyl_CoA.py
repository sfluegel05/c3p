"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
from rdkit import Chem

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Coenzyme A structure patterns using SMARTS
    # This pattern detects phosphate groups and ribose sugar components
    coenzyme_a_pattern = Chem.MolFromSmarts("O=P(OC[C@@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12)O")
    if not mol.HasSubstructMatch(coenzyme_a_pattern):
        return False, "Coenzyme A structure not found"
    
    # Thioester linkage pattern (R-C(=O)-S-R)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage not found"
    
    # Verify presence of a fatty acid hydrocarbon chain
    # Fatty acids typically have long hydrocarbon chains (e.g., C12+)
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("CC[CH2]CCCCCCCCCC")
    if not mol.HasSubstructMatch(hydrocarbon_chain_pattern):
        return False, "Long fatty acid hydrocarbon chain not found"

    # If all patterns are found, it is a fatty acyl-CoA
    return True, "Contains Coenzyme A structure, thioester linkage, and fatty acid moiety"