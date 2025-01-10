"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA derivative based on its SMILES string.

    A 2-enoyl-CoA typically features a CoA moiety, recognized by an adenosine-linked phosphate chain,
    connected by a thioester linkage to an acyl chain with an enoyl group (C=C).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # CoA moiety pattern (focusing on key elements like phosphates and adenosine)
    coa_pattern = Chem.MolFromSmarts("O=P(O)(O)OC[C@H]1O[C@H](COP(=O)(O)O[C@H]2[C@H]([C@H](O)[C@@H]([C@H]2OP(=O)(O)O)n3cnc4c(N)ncnc34)O)C(=O)C=C")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Enoyl group pattern (C=C double bond in acyl chain)
    enoyl_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(enoyl_pattern):
        return False, "No enoyl group found"

    # Thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    return True, "Contains CoA moiety, enoyl group, and thioester linkage"