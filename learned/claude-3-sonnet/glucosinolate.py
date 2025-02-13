"""
Classifies: CHEBI:24279 glucosinolate
"""
"""
Classifies: glucosinolate compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a glucosinolate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the core glucosinolate structure with variations:
    # Central carbon with S and N, allowing for charged and uncharged sulfate
    core_pattern = Chem.MolFromSmarts('[#16]-[#6](=,:[#7]-[#8]-[#16](=[#8])(=[#8])[#8,#8-])-[#6,#7,#8]')
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing core glucosinolate structure"

    # Check for pyranose sugar (glucose) pattern with correct connectivity
    # More specific pattern for β-D-glucopyranose
    sugar_pattern = Chem.MolFromSmarts('[#6]1-[#8]-[#6](-[#16])[#6](-[#8])[#6](-[#8])[#6](-[#8])[#6]1-[#6]-[#8]')
    if not mol.GetSubstructMatches(sugar_pattern):
        return False, "Missing or incorrect glucose moiety"

    # Verify the presence of sulfate group connected to oxime nitrogen
    # Allow for both charged and uncharged forms
    sulfate_pattern = Chem.MolFromSmarts('[#7]-[#8]-[#16](=[#8])(=[#8])[#8,#8-]')
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "Missing or incorrect sulfate group"

    # Count key atoms to ensure reasonable composition
    s_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[#16]')))
    o_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[#8]')))
    n_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[#7]')))
    
    if s_count < 2:  # Need at least 2 sulfur atoms (thioglycoside and sulfate)
        return False, "Insufficient sulfur atoms"
    if o_count < 7:  # Need multiple oxygens for sugar + sulfate
        return False, "Insufficient oxygen atoms"
    if n_count != 1:  # Need exactly 1 nitrogen
        return False, "Incorrect number of nitrogen atoms"

    # Check for hydroxyl groups on sugar
    sugar_oh_pattern = Chem.MolFromSmarts('[#6]1-[#8]-[#6]-[#6](-[#8][H])[#6](-[#8][H])[#6](-[#8][H])[#6]1-[#6]-[#8][H]')
    if not mol.HasSubstructMatch(sugar_oh_pattern):
        return False, "Missing required hydroxyl groups on glucose"

    # Verify C=N bond presence
    cn_double_bond = Chem.MolFromSmarts('[#6]=[#7]')
    if not mol.HasSubstructMatch(cn_double_bond):
        return False, "Missing C=N double bond"

    return True, "Contains complete glucosinolate structure with correct connectivity"