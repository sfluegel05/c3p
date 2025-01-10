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

    # Verify existence of a phosphatidylglycerol structure
    glycerol_pattern = Chem.MolFromSmarts("C(COP(O)(O)=O)OC(=O)C")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No recognizable phosphatidylglycerol backbone found"
    
    # Inositol ring pattern
    # Look for a six-membered ring with hydroxyl groups (R groups can be oxygens for phosphates)
    inositol_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "Inositol ring not found"

    # Phosphate groups directly attached to the inositol ring
    phosphate_pattern = Chem.MolFromSmarts("C(OP(=O)(O)O)C")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    
    # Ensure at least one phosphate group is found
    if len(phosphate_matches) == 0:
        return False, "No phosphorylated inositol ring found"
    
    return True, "Molecule has a phosphatidylinositol backbone with one or more phosphate groups on the inositol ring"