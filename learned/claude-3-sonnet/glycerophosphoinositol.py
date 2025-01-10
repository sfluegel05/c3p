"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (True if molecule is a glycerophosphoinositol, reason for classification)
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("[OX2][P](=[O])([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for inositol group (cyclohexane with 6 OH groups)
    inositol_pattern = Chem.MolFromSmarts("[C@H]1([OH])[C@H]([OH])[C@H]([OH])[C@H]([OH])[C@H]([OH])[C@H]1[OH]")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol group found"

    # Check for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for connection between phosphate and inositol
    phosphoinositol_pattern = Chem.MolFromSmarts("[C][OH][P](=[O])([OH])[O][C]1[C][C][C][C][C]1")
    if not mol.HasSubstructMatch(phosphoinositol_pattern):
        return False, "Phosphate not connected to inositol"

    # Check for at least one ester group (fatty acid attachment)
    ester_pattern = Chem.MolFromSmarts("[#6][CX3](=[OX1])[OX2][#6]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester groups found"

    # Count oxygen atoms (should be at least 11: 6 from inositol, 3 from phosphate, 2+ from esters)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 11:
        return False, f"Too few oxygen atoms ({o_count}) for glycerophosphoinositol"

    # Count carbon atoms (should be at least 9: 6 from inositol, 3 from glycerol)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 9:
        return False, f"Too few carbon atoms ({c_count}) for glycerophosphoinositol"

    # Check for phosphate-glycerol connection
    phosphoglycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4][OX2][P]")
    if not mol.HasSubstructMatch(phosphoglycerol_pattern):
        return False, "Phosphate not connected to glycerol backbone"

    return True, "Contains glycerol backbone with phosphate-linked inositol and fatty acid chains"