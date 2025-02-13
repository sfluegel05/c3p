"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
"""
Classifies: phosphatidylinositol compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol has:
    - An inositol ring (cyclohexane with 6 OH groups)
    - Connected via phosphate to a glycerol backbone
    - Two fatty acid chains attached to the glycerol
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a phosphatidylinositol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for inositol ring pattern - more flexible pattern to match different representations
    # Cyclohexane ring with 6 OH groups in any stereochemistry
    inositol_pattern = Chem.MolFromSmarts("[C@@,C@,CH1,CH2]1[C@@,C@,CH1,CH2]([OH1,OH0])[C@@,C@,CH1,CH2]([OH1,OH0])[C@@,C@,CH1,CH2]([OH1,OH0])[C@@,C@,CH1,CH2]([OH1,OH0])[C@@,C@,CH1,CH2]([OH1,OH0])O1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"
    
    # Look for phosphate group connected to inositol
    # More general pattern that captures different representations
    phosphate_pattern = Chem.MolFromSmarts("O[P](=O)(O)OC")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Look for glycerol backbone - more flexible pattern
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4,CH2][CH1X4,CH1][CH2X4,CH2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for two ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found only {len(ester_matches)} ester groups, need at least 2"
    
    # Check for fatty acid chains - more comprehensive pattern
    # Matches both saturated and unsaturated chains
    fatty_acid_pattern = Chem.MolFromSmarts("([CH2][CH2][CH2][CH2])|([CH2][CH][CH][CH2])")
    chain_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(chain_matches) < 2:
        return False, "Missing long carbon chains (fatty acids)"
    
    # Verify phosphate-inositol-glycerol connectivity
    # More flexible pattern that captures different representations
    connectivity_pattern = Chem.MolFromSmarts("O1[C@@,C@,CH1,CH2]([C@@,C@,CH1,CH2]([OH1,OH0])[C@@,C@,CH1,CH2]([OH1,OH0])[C@@,C@,CH1,CH2]([OH1,OH0])[C@@,C@,CH1,CH2]([OH1,OH0])[C@@,C@,CH1,CH2]1[OH1,OH0])OP(=O)(O)OC")
    if not mol.HasSubstructMatch(connectivity_pattern):
        return False, "Phosphate not correctly connected to inositol and glycerol"
    
    # Count key atoms to verify overall composition
    phosphorus_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if phosphorus_count != 1:
        return False, f"Should have exactly 1 phosphorus atom, found {phosphorus_count}"
    
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 12:  # minimum expected: 6 OH on inositol + 3 phosphate + 2 ester + 1 glycerol
        return False, f"Too few oxygen atoms for phosphatidylinositol (found {oxygen_count}, need at least 12)"
    
    return True, "Contains inositol ring connected via phosphate to glycerol backbone with two fatty acid chains"