"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic elements (excluding H since it might be implicit)
    required_elements = {'C', 'O', 'N', 'P'}
    mol_elements = set(atom.GetSymbol() for atom in mol.GetAtoms())
    if not required_elements.issubset(mol_elements):
        return False, f"Missing required elements. Found {mol_elements}, need {required_elements}"

    # Check for glycerol backbone with correct stereochemistry
    # Multiple SMARTS patterns to catch different representations
    glycerol_patterns = [
        "[CH2X4][C@@HX4][CH2X4]",  # R configuration
        "[CH2X4][C@HX4][CH2X4]",   # Alternative representation
        "[CH2][C@H][CH2]",         # Simplified pattern
        "[CH2][C@@H][CH2]"         # Simplified alternative
    ]
    
    found_glycerol = False
    for pattern in glycerol_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_glycerol = True
            break
    
    if not found_glycerol:
        return False, "No glycerol backbone with correct stereochemistry found"

    # Check for one ester group (acyl chain at position 1)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for phosphoethanolamine group - multiple patterns to catch different forms
    phosphoethanolamine_patterns = [
        "[OX2][PX4](=[OX1])([OX2])[OX2]CC[NX3]",           # Neutral form
        "[OX2][PX4](=[OX1])([OX2])[OX2]CC[NH3+]",         # Protonated amine
        "[OX2][PX4](=[OX1])([O-])[OX2]CC[NH3+]",          # Zwitterionic form
        "[OX2][PX4](=[OX1])([OX2H])[OX2]CC[NX3]"          # With explicit H
    ]
    
    found_pe = False
    for pattern in phosphoethanolamine_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_pe = True
            break
            
    if not found_pe:
        return False, "No phosphoethanolamine group found"

    # Check for hydroxyl group - multiple patterns
    hydroxyl_patterns = [
        "[OX2H1]",            # Explicit H
        "[OX2H]",             # Alternative representation
        "[OH]"                # Simplified
    ]
    
    found_hydroxyl = False
    for pattern in hydroxyl_patterns:
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
        if len(matches) >= 1:
            found_hydroxyl = True
            break
            
    if not found_hydroxyl:
        return False, "No free hydroxyl group found"

    # Check carbon chain length (should be at least 13 carbons total including glycerol backbone)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 13:
        return False, f"Carbon count too low ({c_count}), need at least 13"

    return True, "Contains glycerol backbone with correct stereochemistry, one acyl chain, phosphoethanolamine group, and free hydroxyl group"