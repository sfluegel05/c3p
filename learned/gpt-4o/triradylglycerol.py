"""
Classifies: CHEBI:76579 triradylglycerol
"""
from rdkit import Chem

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triradylglycerol, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Pattern to identify glycerol backbone with three substituents
    glycerol_pattern = Chem.MolFromSmarts("C(O)C(O)C(O)")
    matching_atoms = mol.GetSubstructMatches(glycerol_pattern)
    if not matching_atoms:
        return False, "No glycerol backbone found with correct substituents"
    
    # Pattern to identify ester and ether bonds
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ether_pattern = Chem.MolFromSmarts("COC")
    
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    
    # Check for at least three possible linkage groups (ester or ether)
    if len(ester_matches) + len(ether_matches) < 3:
        return False, f"Found {len(ester_matches) + len(ether_matches)} ester or ether linkages, need at least 3"

    acyl_groups = set()
    for match in ester_matches:
        acyl_groups.add(mol.GetBondWithIdx(match[1]).GetBeginAtom().GetSymbol())
    
    # Test for diversity in substituents
    enough_diversity = False
    if len(acyl_groups) >= 2:
        enough_diversity = True
    
    # Check if substituents have enough diversity and presence
    if enough_diversity or len(ether_matches) > 1:
        return True, "Is a triradylglycerol with enough diversity in substituents"
    
    return False, "Insufficient diversity among substituents"