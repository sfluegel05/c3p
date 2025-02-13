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

    # Check for glycerol backbone with specific substitution pattern
    # [CH2]-[CH]-[CH2] backbone
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for phosphoserine group
    # Phosphoserine moiety connected to glycerol
    phosphoserine_pattern = Chem.MolFromSmarts("[PX4](=O)([OH])([OH])[OX2][CH2][CH]([NH2])[CX3](=O)[OH]")
    if not mol.HasSubstructMatches(phosphoserine_pattern):
        return False, "No phosphoserine group found"

    # Check for acyl group at sn-1 position
    # Look for -O-C(=O)-[C,H] pattern that's not part of the serine carboxyl
    acyl_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6;!$([CX3](=O)[OH])]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 1:
        return False, f"Found {len(acyl_matches)} acyl groups, need exactly 1"

    # Check for middle carbon with hydroxyl group (sn-2 position)
    sn2_hydroxyl = Chem.MolFromSmarts("[CH2X4][CH]([OH])[CH2X4]")
    if not mol.HasSubstructMatch(sn2_hydroxyl):
        return False, "No free hydroxyl group at sn-2 position"

    # Check for fatty acid chain (at least 4 carbons in chain)
    fatty_chain = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(fatty_chain):
        return False, "No fatty acid chain found"

    # Count key atoms to verify overall composition
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if p_count != 1:
        return False, "Must have exactly one phosphorus atom"
    if n_count != 1:
        return False, "Must have exactly one nitrogen atom"
    if o_count < 8:  # minimum number of oxygens needed
        return False, "Insufficient number of oxygen atoms"

    # Verify complete structure connectivity
    # The phosphoserine should be connected to glycerol backbone
    complete_pattern = Chem.MolFromSmarts("[CH2X4][CH]([OH])[CH2X4][OX2][PX4](=O)([OH])[OX2][CH2][CH]([NH2])[CX3](=O)[OH]")
    if not mol.HasSubstructMatch(complete_pattern):
        return False, "Incorrect connectivity of phosphoserine to glycerol backbone"

    return True, "Contains glycerol backbone with single acyl chain at sn-1 position and phosphoserine group at sn-3 position"