"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
from rdkit import Chem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    This involves checking for a CoA moiety, 3-oxo group, and a thioester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-oxo-fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General pattern for the nucleotide part of CoA, considering flexibility in phosphate linkage
    coa_pattern = Chem.MolFromSmarts("COP(=O)(O)[O;D1]C[C@@H]1O[C@H]([C@H](O)[C@@H]1O)n2cnc3[nH]c[nH]c23")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"
    
    # Broader pattern for the 3-oxo group, now allowing more variety in adjacent groups
    oxo_group_pattern = Chem.MolFromSmarts("C(CC(=O))[!O]=O")
    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "No 3-oxo-fatty acid group found"

    # Thioester linkage 'C(=O)S' often observed in fatty acyl-CoA conjugates
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    return True, "Structure matches 3-oxo-fatty acyl-CoA with core CoA moiety, 3-oxo group, and thioester linkage"