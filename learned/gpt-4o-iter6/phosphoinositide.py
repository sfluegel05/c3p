"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol with one or more phosphate groups on the inositol ring.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        (bool, str): Tuple where first element indicates if the molecule is a phosphoinositide,
                      and the second element provides the reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone with fatty acids (ester connections)
    glycerol_pattern = Chem.MolFromSmarts("[CX3](=O)O[CH2][CX2]([O])")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with fatty acids found"

    # Inositol ring detection: Allow for variations in SMILES, capture stereochemistry flexibility
    inositol_pattern = Chem.MolFromSmarts("C1([C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O1)O)")

    # If no match, try without stereochemistry, as some databases use simplified notations
    if not mol.HasSubstructMatch(inositol_pattern):
        inositol_pattern_non_chiral = Chem.MolFromSmarts("C1(C(O)C(O)C(O)C(O)C1O)")
        if not mol.HasSubstructMatch(inositol_pattern_non_chiral):
            return False, "Inositol ring not found"

    # Check for phosphate groups bonded to the inositol ring (at least one required)
    phosphate_pattern = Chem.MolFromSmarts("C(OP(=O)(O)O)C")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphorylated inositol ring found"

    return True, "Molecule has a phosphatidylinositol backbone with one or more phosphate groups on the inositol ring"