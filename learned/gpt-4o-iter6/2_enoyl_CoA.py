"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA is defined as an unsaturated fatty acyl-CoA in which the S-acyl group contains a double bond between positions 2 and 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # More comprehensive pattern capturing the C(=O)S-C=C motif, allowing flexibility for substitution and stereochemistry
    enoyl_pattern = Chem.MolFromSmarts("C(=O)SC([#6&v4])=C([#6])")  # Flexibility in stereochemistry and substituents
    
    # Comprehensive pattern for CoA based on characteristic phosphopantetheine structure
    coA_pattern = Chem.MolFromSmarts("COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O[C@@H]1O)O)n2cnc3c(N)ncnc23")

    # Check for the 2-enoyl pattern
    if not mol.HasSubstructMatch(enoyl_pattern):
        return False, "No 2-enoyl pattern (enoyl thioester linkage with correct positioning) found"

    # Ensure the Coenzyme A (CoA) structure is present
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "No Coenzyme A backbone detected"

    return True, "Contains 2-enoyl feature with Coenzyme A structure"