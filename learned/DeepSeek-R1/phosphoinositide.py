"""
Classifies: CHEBI:18179 phosphoinositide
"""
"""
Classifies: CHEBI_xxxxx phosphoinositide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol with at least one additional phosphate on the inositol ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for phosphatidylinositol core structure: glycerol with two esters and a phosphate connected to inositol
    # Glycerol pattern: C-C-C with two ester groups and one phosphate group
    # Inositol is a cyclohexane with multiple hydroxyls connected via phosphate
    pi_pattern = Chem.MolFromSmarts(
        "[CH2]-[CH](-[OX2]-C(=O)-[!O])-[CH2]-[OX2]-P(=O)(-[OX2])-[OX2][C]1(-[OH])[C](-[OH])[C](-[OH])[C](-[OH])[C](-[OH])[C]1(-[OH])"
    )
    if not mol.HasSubstructMatch(pi_pattern):
        return False, "Not a phosphatidylinositol"
    
    # Check for additional phosphate groups on the inositol ring
    # Look for O-P groups attached to the inositol (excluding the connecting phosphate)
    inositol_phosphate_pattern = Chem.MolFromSmarts("[C](-[OH])-O-P(=O)([OX2])-[OX2]")
    if mol.HasSubstructMatch(inositol_phosphate_pattern):
        return True, "Phosphorylated inositol detected"
    
    return False, "No additional phosphates on inositol"