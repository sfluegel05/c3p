"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
    based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches the structure, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for stereochemistry and backbone pattern with the glycerol core structure
    glycerol_sn_pattern = Chem.MolFromSmarts("[C@H]([CH2]O[C@@H])(O)CO")
    if not mol.HasSubstructMatch(glycerol_sn_pattern):
        return False, "No correct glycerol (sn) backbone found"
    
    # Check for ether linkage indicating an alkyl chain at position 1
    alkyl_ether_pattern = Chem.MolFromSmarts("OCC")  # Adjusted for alkyl O linkage
    if not mol.HasSubstructMatch(alkyl_ether_pattern):
        return False, "No ether linkage indicating alkyl chain found"

    # Check for ester linkage indicating an acyl chain at position 2
    ester_acyl_pattern = Chem.MolFromSmarts("OC(=O)C")  # Acyl pattern checker
    if not mol.HasSubstructMatch(ester_acyl_pattern):
        return False, "No ester linkage indicating acyl chain found"
    
    # Verify presence of phosphocholine group with appropriate linkage
    phosphocholine_pattern = Chem.MolFromSmarts("N(C)(C)CCOP(=O)([O-])O")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No complete phosphocholine moiety found"

    # If all matches are found, it is a correct structure
    return True, "Structure matches 2-acyl-1-alkyl-sn-glycero-3-phosphocholine"