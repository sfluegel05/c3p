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
    Must have an inositol ring with at least one additional phosphate group beyond
    the one connecting to glycerol/lipid.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for inositol ring core (cyclohexane with 6 oxygens)
    inositol_pattern = Chem.MolFromSmarts("[C]1([O])[C]([O])[C]([O])[C]([O])[C]([O])[C]1[O]")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Find all phosphate groups attached to the inositol ring
    # This pattern matches phosphates directly on the ring
    inositol_phosphate = Chem.MolFromSmarts("[C]1[C][C][C][C][C]1OP(=O)([O,OH])[O,OH]")
    phosphate_matches = mol.GetSubstructMatches(inositol_phosphate)
    
    if not phosphate_matches:
        return False, "No phosphate groups attached to inositol"

    # Count phosphates that are NOT connecting to a carbon chain (i.e., not the glycerol-connecting phosphate)
    additional_phosphate = Chem.MolFromSmarts("[C]1[C][C][C][C][C]1OP(=O)([O,OH])[O,OH]")
    phosphate_count = len(mol.GetSubstructMatches(additional_phosphate))
    
    # Look for glycerol/lipid attachment via phosphate
    glycerol_phosphate = Chem.MolFromSmarts("[CH2X4]-[CHX4]-[CH2X4]-OP(=O)([O,OH])-O")
    lipid_phosphate = Chem.MolFromSmarts("COP(=O)([O,OH])O[C]1[C][C][C][C][C]1")
    
    has_glycerol = mol.HasSubstructMatch(glycerol_phosphate)
    has_lipid = mol.HasSubstructMatch(lipid_phosphate)
    
    if not (has_glycerol or has_lipid):
        # Special case: Check for modified phosphoinositides without typical glycerol
        modified_phosphate_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C][C]1OP(=O)(O)O")
        if mol.GetSubstructMatches(modified_phosphate_pattern):
            if phosphate_count >= 2:  # Need at least 2 phosphates for modified structures
                return True, f"Modified phosphoinositide with {phosphate_count} phosphate groups"
        return False, "No phospholipid attachment found"

    # For standard phosphoinositides, we need at least one additional phosphate
    # beyond the one connecting to glycerol/lipid
    if phosphate_count <= 1:
        return False, "Only has connecting phosphate, no additional phosphorylation"

    # Look for ester groups (fatty acid attachments) for standard phosphoinositides
    if has_glycerol:
        ester_pattern = Chem.MolFromSmarts("C(=O)OC")
        ester_matches = mol.GetSubstructMatches(ester_pattern)
        if len(ester_matches) < 1:
            return False, "No ester groups found"

    return True, f"Phosphoinositide with {phosphate_count} phosphate group(s) on inositol ring"