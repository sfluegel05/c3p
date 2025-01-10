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
    Must have:
    - Glycerol backbone with R configuration
    - Acyl group at position 1
    - Free hydroxyl at position 2
    - Phosphoethanolamine at position 3
    
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

    # Check for required elements
    required_elements = {'C', 'O', 'N', 'P'}
    mol_elements = set(atom.GetSymbol() for atom in mol.GetAtoms())
    if not required_elements.issubset(mol_elements):
        return False, f"Missing required elements. Found {mol_elements}, need {required_elements}"

    # Basic glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CH1X4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for ester group at position 1
    ester_pattern = Chem.MolFromSmarts("[CH2X4]OC(=O)[#6]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group at position 1"

    # Check for hydroxyl at position 2
    hydroxyl_pattern = Chem.MolFromSmarts("[CH1X4][OX2H1]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No free hydroxyl group at position 2"

    # Check for phosphoethanolamine group - try both neutral and charged forms
    pe_patterns = [
        Chem.MolFromSmarts("COP(=O)(O)OCCN"),  # Neutral form
        Chem.MolFromSmarts("COP(=O)([O-])OCC[NH3+]"),  # Zwitterionic form
    ]
    
    found_pe = False
    for pattern in pe_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            found_pe = True
            break
            
    if not found_pe:
        return False, "No phosphoethanolamine group found"

    # Verify no acylation at position 2
    wrong_acyl = Chem.MolFromSmarts("[CH1X4](OC(=O))[CH2X4]")
    if wrong_acyl is not None and mol.HasSubstructMatch(wrong_acyl):
        return False, "Has acyl group at position 2"

    # Check for reasonable chain length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 7:
        return False, f"Carbon count too low ({c_count}), need at least 7"

    # Check for expected number of phosphorus atoms
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 1:
        return False, f"Expected exactly 1 phosphorus atom, found {p_count}"

    # Verify basic connectivity
    # Should have: CH2-CH-CH2 backbone with:
    # - O-acyl at position 1
    # - OH at position 2
    # - O-P at position 3
    basic_connectivity = Chem.MolFromSmarts("[CH2X4]OC(=O)[#6][CH1X4]([OX2H1])[CH2X4]OP")
    if not mol.HasSubstructMatch(basic_connectivity):
        return False, "Core structural connectivity not found"

    return True, "Valid 1-acyl-sn-glycero-3-phosphoethanolamine structure found"