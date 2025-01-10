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
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for inositol ring (cyclohexane with oxygen substituents)
    # More flexible pattern that matches cyclohexane with O-substitutions
    inositol_pattern = Chem.MolFromSmarts("[C]1([O,OH])[C]([O,OH])[C]([O,OH])[C]([O,OH])[C]([O,OH])[C]1[O,OH]")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Look for phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("O-P(=O)(-[O,OH])-[O,OH]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphate groups found"

    # Look for glycerol backbone connected to phosphate
    glycerol_phosphate_pattern = Chem.MolFromSmarts("[CH2X4]-[CHX4]-[CH2X4]-O-P(=O)(-[O,OH])-O")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No phosphate-linked glycerol backbone found"

    # Look for ester groups (fatty acid attachments)
    ester_pattern = Chem.MolFromSmarts("C(=O)-O-[CH2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, "Missing fatty acid chains"

    # Verify phosphate connection to inositol
    inositol_phosphate = Chem.MolFromSmarts("[C]1([O,OH])[C]([O,OH])[C]([O,OH])[C]([O,OH])[C]([O,OH])[C]1[O]-P(=O)(-[O,OH])-O")
    if not mol.HasSubstructMatch(inositol_phosphate):
        return False, "No phosphate group connected to inositol"

    # Count carbons in longest chain to verify fatty acids
    carbon_chain = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "No long carbon chains found"

    # Additional check for phosphorylation on inositol ring
    inositol_phosphate_direct = Chem.MolFromSmarts("[C]1[C][C][C][C][C]1-O-P(=O)(-[O,OH])-[O,OH]")
    if not mol.HasSubstructMatch(inositol_phosphate_direct):
        return False, "No phosphorylation on inositol ring"

    return True, "Contains phosphorylated inositol ring with glycerol-phosphate backbone and fatty acid chains"