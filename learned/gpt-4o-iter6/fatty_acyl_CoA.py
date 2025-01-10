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

    # Broader Coenzyme A structure including adenosine and phosphate groups
    coenzyme_a_pattern = Chem.MolFromSmarts("O=P(OC[C@@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)O)n1cnc2c(N)ncnc12)O")
    if not mol.HasSubstructMatch(coenzyme_a_pattern):
        return False, "Coenzyme A structure not found"
    
    # Thioester linkage pattern (R-C(=O)-S-CoA)
    thioester_pattern = Chem.MolFromSmarts("C(=O)SSCCNC(=O)C")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage with CoA not found"
    
    # More flexible pattern for a fatty acid hydrocarbon chain to allow for variations
    # Matches a series of carbon atoms which may include some branching and double bonds
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("C-C-C-C-C-C-C-C")  # Representing an ideal fluid portion
    if not mol.HasSubstructMatch(hydrocarbon_chain_pattern):
        return False, "No sufficient hydrocarbon chain part detected"

    # Passed all checks - return true with reason
    return True, "Contains Coenzyme A structure, thioester linkage, and fatty acid moiety"