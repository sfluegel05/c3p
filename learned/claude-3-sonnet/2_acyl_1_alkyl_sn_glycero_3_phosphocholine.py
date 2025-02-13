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

    # Check for complete phosphocholine group structure
    # More specific SMARTS pattern to avoid partial matches
    phosphocholine = Chem.MolFromSmarts("[O][P](=[O])([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine):
        return False, "No phosphocholine group found"
    
    # Count phosphocholine groups - should be exactly one
    phosphocholine_matches = len(mol.GetSubstructMatches(phosphocholine))
    if phosphocholine_matches != 1:
        return False, f"Found {phosphocholine_matches} phosphocholine groups, need exactly 1"

    # Check for the complete glycerol backbone with correct substitution pattern
    # [C@H] ensures the sn-glycero configuration
    # Position 1: -CH2-O-R (ether)
    # Position 2: -CH-O-C(=O)- (ester)
    # Position 3: -CH2-O-P (phosphate)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4]-[OX2]-[CX4H2,CX4H3].[CH2X4]-[OX2]-[PX4].[C@H](-[OX2]-[CX3](=O))")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Missing correct glycerol backbone structure with required substitution pattern"

    # Check for ether linkage at position 1 (long chain)
    ether_pattern = Chem.MolFromSmarts("[CX4H2]-[OX2]-[CX4H2]-[CX4H2]-[CX4H2]")
    if not mol.HasSubstructMatch(ether_pattern):
        return False, "No long-chain alkyl ether found at position 1"

    # Check for ester group at position 2
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])-[CX4]")
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    if ester_matches != 1:
        return False, f"Found {ester_matches} ester groups, need exactly 1"

    # Verify minimum size requirements
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 15:  # Minimum size considering the basic structure
        return False, "Molecule too small for expected structure"

    # Check for required number of oxygens (ether, ester, phosphate)
    n_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if n_oxygens < 7:  # 1 ether + 2 ester + 4 phosphate
        return False, "Insufficient oxygen atoms for required functional groups"

    return True, "Contains glycerol backbone with alkyl ether at position 1, acyl group at position 2, and phosphocholine at position 3"