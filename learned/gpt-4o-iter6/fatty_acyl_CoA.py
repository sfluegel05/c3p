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

    # Broader Coenzyme A structure including adenosine, thiol, and phosphate groups
    coenzyme_a_pattern = Chem.MolFromSmarts("O=S(=O)(O[C@H]1[C@H]([C@@H]([C@H](O1)n2cnc3c(ncnc23)N)OP(=O)(O)O)OP(=O)(O)O)")
    if not mol.HasSubstructMatch(coenzyme_a_pattern):
        return False, "Coenzyme A structure not found"
    
    # Revised thioester linkage pattern (R-C(=O)-S-CoA)
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)C")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage with CoA not found"
    
    # Allowing for double bonds in a flexible fatty acid hydrocarbon chain
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]")
    if not mol.HasSubstructMatch(hydrocarbon_chain_pattern):
        return False, "No sufficient hydrocarbon chain part detected"

    # Passed all checks - return true with reason
    return True, "Contains Coenzyme A structure, thioester linkage, and fatty acid moiety"