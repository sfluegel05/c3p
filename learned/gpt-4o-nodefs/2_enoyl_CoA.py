"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA derivative based on its SMILES string.

    A 2-enoyl-CoA derivative features a CoA moiety (Coenzyme A), typically recognized by
    an adenosine-linked phosphate chain, linked via thioester to a fatty acyl chain with an enoyl group (C=C).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) - True if molecule is a 2-enoyl-CoA derivative with a reason, False otherwise
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Simplified CoA part pattern
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)O")
    if not coa_pattern or not mol.HasSubstructMatch(coa_pattern):
        return False, "No or incomplete CoA moiety found"
    
    # Enoyl group pattern (C=C double bond in acyl chain)
    enoyl_pattern = Chem.MolFromSmarts("C=C")
    if not enoyl_pattern or not mol.HasSubstructMatch(enoyl_pattern):
        return False, "No enoyl group found"

    # Thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not thioester_pattern or not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Pass all checks
    return True, "Contains CoA moiety, enoyl group, and thioester linkage"