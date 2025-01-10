"""
Classifies: CHEBI:17389 2-monoglyceride
"""
from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride has a glycerol backbone with an acyl substituent at the 2-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to improve substructure matching
    mol = Chem.AddHs(mol)

    # Look for the specific glycerol pattern where the middle hydroxyl group is esterified
    # Glycerol: HO-CH2-CH(OH)-CH2OH, 2-monoacyl: HO-CH2-CH(O-C(=O)-R)-CH2OH
    glycerol_pattern = Chem.MolFromSmarts("O[C@@H](CO)C(=O)[*]")
    if glycerol_pattern is None:
        return False, "Error in constructing SMARTS pattern"

    # Check for the glycerol backbone with esterified acyl at the correct position
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No acyl substituent at the 2-position via ester linkage"

    # Exclude entities with non-trivial additional substituents characteristic of complex molecules
    complex_substituents = Chem.MolFromSmarts("[#7,#15,#16]")  # Exclude nitrogen, phosphorus, sulfur
    if mol.HasSubstructMatch(complex_substituents):
        return False, "Contains complex substituents that are not characteristic of a 2-monoglyceride"

    return True, "Contains glycerol backbone with an acyl substituent at the 2-position (ester linkage)"