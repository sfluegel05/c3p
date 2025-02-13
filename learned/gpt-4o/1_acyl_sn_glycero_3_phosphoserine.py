"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    This includes a glycerol backbone with an acyl chain at the 1-hydroxy position and a phosphoserine group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoserine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone and phosphoserine structure with correct stereochemistry
    backbone_pattern = Chem.MolFromSmarts("O[C@@H](CO[P](=O)(O)OC[C@H](N)C(=O)O)CO")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No glycerol backbone with correct stereochemistry found"
    
    # Ensure the presence of a single acyl chain with ester linkage where expected
    acyl_pattern = Chem.MolFromSmarts("C(=O)OC[C@H](O)CO[P]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 1:
        return False, f"Expected exactly 1 acyl chain bound to the glycerol 1-hydroxy position, found {len(acyl_matches)}"
    
    # Ensure there are no additional acyl chains (avoiding two-chain phospholipids)
    # Check for ester linkage and count them
    ester_link_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_link_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester linkages, expected exactly 1"

    return True, "Contains glycerol backbone with acyl chain and attached phosphoserine group"