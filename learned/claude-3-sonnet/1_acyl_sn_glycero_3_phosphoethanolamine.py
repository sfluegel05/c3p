"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

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

    # Check for glycerol backbone with correct stereochemistry
    # [C@@H] indicates R configuration at position 2
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][C@@HX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        # Try alternative SMARTS for different representation of same stereochemistry
        glycerol_pattern_alt = Chem.MolFromSmarts("[CH2X4][C@HX4][CH2X4]")
        if not mol.HasSubstructMatch(glycerol_pattern_alt):
            return False, "No glycerol backbone with correct stereochemistry found"

    # Check for one ester group (acyl chain at position 1)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for phosphoethanolamine group
    # Pattern matches -O-P(=O)(O)-O-CH2-CH2-N
    phosphoethanolamine = Chem.MolFromSmarts("[OX2][PX4](=[OX1])([OX2])[OX2]CC[NX3]")
    if not mol.HasSubstructMatch(phosphoethanolamine):
        return False, "No phosphoethanolamine group found"

    # Check for hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "No free hydroxyl group found"

    # Verify elements present
    formula = CalcMolFormula(mol)
    required_elements = {'C', 'H', 'O', 'N', 'P'}
    mol_elements = set(atom.GetSymbol() for atom in mol.GetAtoms())
    if not required_elements.issubset(mol_elements):
        return False, f"Missing required elements. Found {mol_elements}, need {required_elements}"

    # Check carbon chain length (should be at least 13 carbons total including glycerol backbone)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 13:
        return False, f"Carbon count too low ({c_count}), need at least 13"

    return True, "Contains glycerol backbone with correct stereochemistry, one acyl chain, phosphoethanolamine group, and free hydroxyl group"