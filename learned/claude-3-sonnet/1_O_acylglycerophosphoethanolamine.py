"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
"""
Classifies: 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone with phosphoethanolamine
    # [CH2X4]-[CHX4]-[CH2X4] where one CH2 has ester, one has phosphate, and CH has OH
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for phosphoethanolamine group
    # -P(=O)(O)-O-CH2-CH2-N
    phosphoethanolamine = Chem.MolFromSmarts("[PX4](=O)([OH,O-])[OX2][CH2X4][CH2X4][NX3,NX4]")
    if not mol.HasSubstructMatch(phosphoethanolamine):
        return False, "No phosphoethanolamine group found"

    # Check for ester group (acyl)
    # R-C(=O)-O-
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester (acyl) group found"
    
    # Count ester groups - should only be one
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    if ester_matches != 1:
        return False, f"Found {ester_matches} ester groups, should be exactly 1"

    # Check for free hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No free hydroxyl group found"

    # Verify basic composition
    # Should have exactly:
    # - One phosphorus
    # - At least 5 oxygens (1 carbonyl, 1 ester, 1 hydroxyl, 2 phosphate)
    # - One nitrogen
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    if p_count != 1:
        return False, f"Should have exactly 1 phosphorus, found {p_count}"
    if o_count < 5:
        return False, f"Should have at least 5 oxygens, found {o_count}"
    if n_count != 1:
        return False, f"Should have exactly 1 nitrogen, found {n_count}"

    # Check for acyl chain length - should be at least 4 carbons
    carbon_chain = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "Acyl chain too short"

    return True, "Contains glycerol backbone with phosphoethanolamine group and single acyl substituent"