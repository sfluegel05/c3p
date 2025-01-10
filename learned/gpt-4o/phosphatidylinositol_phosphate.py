"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
from rdkit import Chem

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Identifies phosphatidylinositol phosphates with a glycerol backbone, an inositol ring,
    and one or more phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: Identification result, True if molecule is a phosphatidylinositol phosphate
        str: Explanation for the classification
    """
    
    # Create rdkit Molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Adjusted glycerol pattern, more generalized
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](COC=O)OC=O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found with two ester bonds"

    # Improved inositol pattern: oxygen and correct stereochemistry (no stereo matching for simplicity)
    inositol_pattern = Chem.MolFromSmarts("C1C(O)C(O)C(O)C(O)C1(O)")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Improved phosphate group pattern: matches both charged and neutral phosphate ester
    phosphate_pattern = Chem.MolFromSmarts("O=P(O)(O)O")
    phosphate_count = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_count == 0:
        return False, "No phosphate groups found"

    return True, "Contains glycerol backbone with inositol ring and one or more phosphate groups"

# Test example
result, reason = is_phosphatidylinositol_phosphate("CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@H]1O)OC(=O)CCCCCCC\C=C/CCCCCCCC")
print(result, reason)  # Expected output: True, with detailed reason