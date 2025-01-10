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

    # Look for inositol ring (cyclohexane with multiple OH groups)
    inositol_pattern = Chem.MolFromSmarts("[OH1]-[CH1]-1-[CH1](-[OH1])-[CH1](-[OH1,OPO3H2])-[CH1](-[OH1,OPO3H2])-[CH1](-[OH1,OPO3H2])-[CH1](-[OH1,OPO3H2])-1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Look for phosphate group connecting glycerol to inositol
    glycerol_phosphate_pattern = Chem.MolFromSmarts("[CH2X4]-[CHX4]-[CH2X4]-O-P(=O)(-O)-O")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No phosphate-linked glycerol backbone found"

    # Look for additional phosphate groups on inositol
    phosphate_pattern = Chem.MolFromSmarts("O-P(=O)(-O)-O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:  # At least 2 phosphates (one connecting glycerol, one on inositol)
        return False, "Insufficient phosphate groups"

    # Look for ester groups (fatty acid attachment points)
    ester_pattern = Chem.MolFromSmarts("[#6]-C(=O)-O-[CH2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, "Missing fatty acid chains"

    # Count phosphorus atoms to verify phosphorylation
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count < 2:  # Must have at least 2 phosphorus atoms
        return False, "Insufficient phosphorus atoms"

    # Verify the presence of long carbon chains (fatty acids)
    carbon_chain = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "No long carbon chains found"

    return True, "Contains inositol ring with multiple phosphate groups and fatty acid chains"