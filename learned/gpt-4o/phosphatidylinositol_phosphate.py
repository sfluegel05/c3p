"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
from rdkit import Chem

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Identifies if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    The essential structural elements include a glycerol backbone, a myo-inositol
    head group, and one or more phosphate groups attached.

    Args:
        smiles (str): SMILES string of the chemical entity

    Returns:
        bool: True if SMILES represents a phosphatidylinositol phosphate, False otherwise
        str: Explanation for classification
    """
    
    # Convert SMILES to RDKit Molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a phosphoglycerol moiety
    phosphoglycerol_pattern = Chem.MolFromSmarts("O[C@H](COP(O)(=O)(O))COC(=O)")
    if not mol.HasSubstructMatch(phosphoglycerol_pattern):
        return False, "No phosphoglycerol subgroup found"

    # Check for a myo-inositol ring
    myo_inositol_pattern = Chem.MolFromSmarts("C1(O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1(O)")
    if not mol.HasSubstructMatch(myo_inositol_pattern):
        return False, "No myo-inositol ring found"

    # Check for phosphate groups, allowing multiple phosphates
    phosphate_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("OP(O)(O)=O"))
    if len(phosphate_groups) < 1:
        return False, "No phosphate groups attached to inositol"

    return True, "Identified as a phosphatidylinositol phosphate with a glycerol backbone, myo-inositol ring, and phosphate groups"

# Test Example
result, reason = is_phosphatidylinositol_phosphate("CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](OP(O)(O)=O)[C@H](O)[C@H]1O)OC(=O)CCCCCCC\C=C/CCCCCCCC")
print(result, reason)