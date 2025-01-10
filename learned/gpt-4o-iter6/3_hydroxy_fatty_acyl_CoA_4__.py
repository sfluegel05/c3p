"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA_4_(smiles: str):
    """
    Classifies a molecule as a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a 3-hydroxy fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Check for a 3-hydroxy group attached to a central carbon, with and without chirality
    hydroxy_chiral_pattern = Chem.MolFromSmarts("[C@H](O)C")
    hydroxy_nonchiral_pattern = Chem.MolFromSmarts("COC")
    
    if not (mol.HasSubstructMatch(hydroxy_chiral_pattern) or mol.HasSubstructMatch(hydroxy_nonchiral_pattern)):
        return False, "No 3-hydroxy group found on fatty acid chain"
    
    # Coenzyme A moiety pattern with thioester linkage
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A moiety with thioester linkage not found"

    # Look for deprotonated phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([O-])[O-]")
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_matches < 2:
        return False, f"Found {phosphate_matches} deprotonated phosphate groups, need at least 2"

    # Fatty acid chain pattern (allowing for some variability in chain length)
    chain_length_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC")  # Minimum pattern ensuring chain length similarity
    if not mol.HasSubstructMatch(chain_length_pattern):
        return False, "Insufficient long fatty acyl chain detected"

    # If all conditions are met
    return True, "The molecule is a 3-hydroxy fatty acyl-CoA(4-)"