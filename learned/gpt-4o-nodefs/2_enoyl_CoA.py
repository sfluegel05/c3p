"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA derivative based on its SMILES string.

    A 2-enoyl-CoA derivative features a specific CoA moiety, linked via thioester to a fatty acyl chain with an enoyl group (C=C).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) - True if molecule is a 2-enoyl-CoA derivative with a reason, False otherwise
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined CoA part pattern
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)[C@@H](C)COP(=O)(O)O[C@H]1[C@@H](O)[C@H](O)C[C@@H]1OP(=O)(O)O[C@@H]2O[C@H](n3cnc4c(N)ncnc34)[C@H](O)[C@H]2O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No complete CoA moiety found"
    
    # Refined enoyl group pattern (more context-specific)
    enoyl_pattern = Chem.MolFromSmarts("C=CC(=O)S")
    if not mol.HasSubstructMatch(enoyl_pattern):
        return False, "No enoyl thioester linkage found"

    # Pass all checks
    return True, "Contains complete CoA moiety and enoyl-thioester linkage"