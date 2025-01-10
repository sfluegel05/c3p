"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: CHEBI:57262 3-sn-phosphatidyl-L-serine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-sn-phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("[OX2][P](=[O])([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for glycerol backbone with specific stereochemistry
    # [C@H] indicates R stereochemistry at carbon 2 of glycerol (sn notation)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][C@H][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with correct stereochemistry found"

    # Check for two ester groups (fatty acid chains)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found only {len(ester_matches)} ester groups, need at least 2"

    # Check for serine moiety
    # Looking for NH2-CH-COOH pattern
    serine_pattern = Chem.MolFromSmarts("[NH2][CH]C(=O)[OH]")
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "No serine moiety found"

    # Count key atoms to verify overall composition
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)  # Phosphorus
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)   # Nitrogen
    
    if p_count != 1:
        return False, "Must have exactly one phosphorus atom"
    if n_count != 1:
        return False, "Must have exactly one nitrogen atom (from serine)"

    # Verify the connection pattern
    # The phosphate should connect glycerol to serine
    phosphoserine_pattern = Chem.MolFromSmarts("[CH2X4]-[OX2]-[P](=[O])([OX2])[OX2]-[CH2][CH]([NH2])C(=O)[OH]")
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "Incorrect connectivity between phosphate and serine"

    # Count carbons in fatty acid chains (should be substantial)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, "Fatty acid chains too short or missing"

    return True, "Contains glycerol backbone with two acyl groups and phosphoserine moiety in correct configuration"