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

    # Pattern for complete 1-acyl-sn-glycero-3-phosphoethanolamine structure
    # [CH2X4] - first carbon with acyl
    # [C@H] - stereocenter with OH
    # [CH2X4] - third carbon with phosphate
    # Includes specific connectivity and phosphoethanolamine group
    core_pattern = """
        [CH2X4;!$(C-[CH2]-[O,N,S])][C@H]([OX2H1])[CH2X4]OP(=[O])(O)OCCN
    """
    
    if not mol.HasSubstructMatch(Chem.MolFromSmarts(core_pattern)):
        return False, "Does not match core 1-acyl-sn-glycero-3-phosphoethanolamine structure"

    # Check specifically for acyl group at position 1
    acyl_pattern = """
        [CH2X4;!$(C-[CH2]-[O,N,S])]OC(=O)[#6]
    """
    if not mol.HasSubstructMatch(Chem.MolFromSmarts(acyl_pattern)):
        return False, "No acyl group at position 1"

    # Check for absence of acylation at position 2
    wrong_acyl_pattern = """
        [C@H](OC(=O)[#6])[CH2X4]
    """
    if mol.HasSubstructMatch(Chem.MolFromSmarts(wrong_acyl_pattern)):
        return False, "Has acyl group at position 2"

    # Check for absence of ether linkages
    ether_pattern = """
        [CH2X4]O[CH2][#6]
    """
    if mol.HasSubstructMatch(Chem.MolFromSmarts(ether_pattern)):
        return False, "Contains ether linkage"

    # Verify phosphoethanolamine group specifically
    pe_patterns = [
        "OP(=O)(O)OCCN",  # Neutral form
        "OP(=O)([O-])OCC[NH3+]",  # Zwitterionic form
    ]
    found_pe = False
    for pattern in pe_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_pe = True
            break
    
    if not found_pe:
        return False, "No valid phosphoethanolamine group found"

    # Count carbons to ensure reasonable chain length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 7:  # Minimum for shortest possible acyl chain
        return False, f"Carbon count too low ({c_count}), need at least 7"

    return True, "Valid 1-acyl-sn-glycero-3-phosphoethanolamine structure found"