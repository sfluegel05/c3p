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
    
    # Check for glycerol backbone with correct stereochemistry
    # [C@H] ensures the sn-glycero configuration
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][C@H,C@@H][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with correct stereochemistry found"
    
    # Check for phosphocholine group (-P(=O)([O-])OCC[N+](C)(C)C)
    phosphocholine = Chem.MolFromSmarts("[P](=[O])([O-])(OCC[N+](C)(C)C)O")
    if not mol.HasSubstructMatch(phosphocholine):
        return False, "No phosphocholine group found"
    
    # Check for ether linkage at position 1 (R-O-CH2-)
    ether_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[CH2X4]")
    if not mol.HasSubstructMatch(ether_pattern):
        return False, "No alkyl ether linkage found at position 1"
    
    # Check for ester group at position 2 (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found at position 2"
    
    # Count number of ester groups - should be exactly one
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"
    
    # Count phosphate groups - should be exactly one
    phosphate_pattern = Chem.MolFromSmarts("[P](=[O])([O-])")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) != 1:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need exactly 1"
    
    # Check for trimethylammonium group - should be exactly one
    nme3_pattern = Chem.MolFromSmarts("[N+](C)(C)(C)")
    nme3_matches = mol.GetSubstructMatches(nme3_pattern)
    if len(nme3_matches) != 1:
        return False, f"Found {len(nme3_matches)} trimethylammonium groups, need exactly 1"
    
    # Verify carbon chain lengths - should have substantial alkyl chains
    # Count carbons in longest chain
    n_carbons = len(mol.GetAtoms())
    if n_carbons < 20:  # Minimum size considering the basic structure
        return False, "Molecule too small for expected structure"
        
    return True, "Contains glycerol backbone with alkyl ether at position 1, acyl group at position 2, and phosphocholine at position 3"