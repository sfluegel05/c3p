"""
Classifies: CHEBI:57643 1,2-diacyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_2_diacyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 1,2-diacyl-sn-glycero-3-phosphocholine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule matches the 1,2-diacyl-sn-glycero-3-phosphocholine structure, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sn-glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("O[C@@H](CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No sn-glycerol backbone found"
        
    # Look for 2 ester groups indicating acyl chains
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"
    
    # Look for phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("P(=O)(O)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "Phosphocholine group not found"

    # Check molecular weight - typical range for such phospholipids is >700 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, "Molecular weight too low for 1,2-diacyl-sn-glycero-3-phosphocholine"

    return True, "Contains sn-glycerol backbone, two acyl chains via ester linkages, and a phosphocholine moiety"