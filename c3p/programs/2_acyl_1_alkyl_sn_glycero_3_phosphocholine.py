"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for complete phosphocholine group with exact connectivity
    phosphocholine = Chem.MolFromSmarts("[O][P](=[O])([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine):
        return False, "No phosphocholine group found"
    
    # Count phosphocholine groups - should be exactly one
    phosphocholine_matches = len(mol.GetSubstructMatches(phosphocholine))
    if phosphocholine_matches != 1:
        return False, f"Found {phosphocholine_matches} phosphocholine groups, need exactly 1"

    # Check for specific glycerol backbone with sn-stereochemistry
    # Note: [C@H] specifies the required stereochemistry for sn-glycero configuration
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4]-[OX2]-[#6][CH2X4]-[C@H](-[OX2]-[CX3]=O)-[CH2X4]-[OX2]-[P]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Missing correct glycerol backbone with sn-stereochemistry"

    # Check for ether linkage specifically at position 1
    ether_pattern = Chem.MolFromSmarts("[CH2X4]-[OX2]-[CX4]")
    if not mol.HasSubstructMatch(ether_pattern):
        return False, "No alkyl ether found at position 1"

    # Check for acyl group specifically at position 2
    acyl_pattern = Chem.MolFromSmarts("[C@H](-[OX2]-[CX3](=[OX1]))-[CH2X4]")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl group found at position 2"

    # Verify phosphocholine at position 3
    position3_pattern = Chem.MolFromSmarts("[C@H](-[CH2X4]-[OX2]-[P])-[OX2]")
    if not mol.HasSubstructMatch(position3_pattern):
        return False, "Incorrect phosphocholine position"

    # Additional checks for complete structure
    # Count carbons to ensure minimum size
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 8:
        return False, "Molecule too small for expected structure"

    # Count oxygens to ensure correct number of functional groups
    n_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if n_oxygens != 6:  # Exactly 6: 1 ether + 2 acyl + 3 phosphate
        return False, "Incorrect number of oxygen atoms"

    # Verify alkyl chain length (minimum 2 carbons)
    alkyl_chain = Chem.MolFromSmarts("[OX2]-[CX4]-[CX4]")
    if not mol.HasSubstructMatch(alkyl_chain):
        return False, "Alkyl chain too short"

    return True, "Contains glycerol backbone with alkyl ether at position 1, acyl group at position 2, and phosphocholine at position 3"