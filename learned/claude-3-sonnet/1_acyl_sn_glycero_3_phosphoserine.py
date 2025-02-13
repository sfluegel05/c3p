"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
"""
Classifies: CHEBI:75944 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoserine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate group connected to serine
    # More flexible pattern that matches various representations
    phosphoserine_pattern = Chem.MolFromSmarts("[PX4](=[OX1])([OX2,OH])([OX2,OH])[OX2][CH2][CH]([NH2,NH3+,N])[CX3](=[OX1])[OX2,OH]")
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "No phosphoserine group found"

    # Check for glycerol backbone 
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.GetSubstructMatches(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for acyl group (ester linkage)
    acyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if not acyl_matches:
        return False, "No acyl group found"

    # Check for fatty acid chain (at least 4 carbons)
    fatty_chain = Chem.MolFromSmarts("C(~C)(~C)~C")
    if not mol.HasSubstructMatch(fatty_chain):
        return False, "No fatty acid chain found"

    # Verify core structure connectivity
    # Phosphate connected to glycerol
    core_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4][OX2][PX4]")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Incorrect connectivity between glycerol and phosphate"

    # Count key atoms to verify overall composition
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if p_count != 1:
        return False, "Must have exactly one phosphorus atom"
    if n_count != 1:
        return False, "Must have exactly one nitrogen atom"
    if o_count < 8:
        return False, "Insufficient number of oxygen atoms"

    # Check for sn-2 hydroxyl group (free OH)
    sn2_pattern = Chem.MolFromSmarts("[CH2X4][CH]([OH,O])[CH2X4]")
    if not mol.HasSubstructMatch(sn2_pattern):
        return False, "Missing hydroxyl group at sn-2 position"

    # Verify the presence of all key components
    if all([
        mol.HasSubstructMatch(phosphoserine_pattern),
        mol.HasSubstructMatch(glycerol_pattern),
        mol.HasSubstructMatch(acyl_pattern),
        mol.HasSubstructMatch(core_pattern)
    ]):
        return True, "Contains glycerol backbone with acyl chain and phosphoserine group in correct positions"

    return False, "Missing one or more required structural features"