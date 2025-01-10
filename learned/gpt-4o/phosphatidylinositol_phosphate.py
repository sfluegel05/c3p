"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    Phosphatidylinositol phosphates are characterized by a glycerol backbone, an inositol ring
    with one or more phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol phosphate, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for glycerol backbone with two ester bonds
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](COC=O)OC=O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found with two ester bonds"

    # Pattern for inositol ring (6-membered oxygen-containing ring with adjacent stereocenters)
    inositol_pattern = Chem.MolFromSmarts("[C@H]1[C@H]([O])C([@H]([O])C1[O])")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Check for phosphate groups (PO4)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)[O-]")
    phosphate_count = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_count < 1:
        return False, "No phosphate groups found"
    
    return True, "Molecule matches phosphatidylinositol phosphate structure with glycerol backbone, inositol ring, and phosphate groups"