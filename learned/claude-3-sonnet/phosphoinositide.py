"""
Classifies: CHEBI:18179 phosphoinositide
"""
"""
Classifies: phosphoinositide
Definition: Any phosphatidylinositol that is phosphorylated at one or more of the hydroxy groups of inositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise
        str: str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for inositol ring core with correct connectivity
    # Match cyclohexane with 6 oxygens (as OH or OR)
    inositol_pattern = Chem.MolFromSmarts("[C]1([O])[C]([O])[C]([O])[C]([O])[C]([O])[C]1[O]")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Look for phosphate groups attached to inositol
    inositol_phosphate = Chem.MolFromSmarts("[C]1([O,OH])[C]([O,OH])[C]([O,OH])[C]([O,OH])[C]([O,OH])[C]1[O]P(=O)([O,OH])[O,OH]")
    phosphate_matches = mol.GetSubstructMatches(inositol_phosphate)
    if not phosphate_matches:
        return False, "No phosphate groups attached to inositol"

    # Look for glycerol-phosphate backbone
    glycerol_phosphate = Chem.MolFromSmarts("[CH2X4]-[CHX4]-[CH2X4]-OP(=O)([O,OH])-O")
    if not mol.HasSubstructMatch(glycerol_phosphate):
        return False, "No glycerol-phosphate backbone found"

    # Count phosphate groups on inositol
    phosphate_on_inositol = Chem.MolFromSmarts("[C]1[C][C][C][C][C]1OP(=O)([O,OH])[O,OH]")
    phosphate_count = len(mol.GetSubstructMatches(phosphate_on_inositol))
    if phosphate_count < 1:
        return False, "No phosphate groups on inositol ring"

    # Look for ester groups (fatty acid attachments)
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester groups found"

    # Verify presence of acyl chains
    acyl_pattern = Chem.MolFromSmarts("CC(=O)O[CH2][CH][CH2]OP")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl chains attached to glycerol backbone"

    # Verify complete phosphoinositide structure:
    # - Inositol ring with at least one phosphate
    # - Connected to glycerol via phosphate
    # - Glycerol with fatty acid chains
    complete_structure = (
        mol.HasSubstructMatch(inositol_phosphate) and
        mol.HasSubstructMatch(glycerol_phosphate) and
        len(ester_matches) >= 1 and
        mol.HasSubstructMatch(acyl_pattern)
    )
    
    if not complete_structure:
        return False, "Incomplete phosphoinositide structure"

    return True, f"Phosphoinositide with {phosphate_count} phosphate group(s) on inositol ring"