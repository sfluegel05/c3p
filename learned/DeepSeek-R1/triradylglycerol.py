"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: CHEBI:17855 triradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.
    A triradylglycerol has a glycerol backbone with three substituents, each being
    acyl (ester), alkyl (ether), or alk-1-enyl (vinyl ether) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triradylglycerol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone with three O-linked substituents
    backbone_pattern = Chem.MolFromSmarts('C(-O)-C(-O)-C(-O)')
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No glycerol backbone with three O-linked groups"

    # Get O atoms from the first backbone match
    matches = mol.GetSubstructMatches(backbone_pattern)
    if not matches:
        return False, "No backbone match found"
    
    # Extract O indices from the first match (C-O positions)
    try:
        o_indices = [matches[0][1], matches[0][3], matches[0][5]]  # C-O positions in SMARTS
    except IndexError:
        return False, "Invalid backbone structure"

    # Define patterns for substituent types
    ester_pattern = Chem.MolFromSmarts('[OX2]-C(=O)')
    vinyl_ether_pattern = Chem.MolFromSmarts('[OX2]-C=C')

    substituent_types = []
    for o_idx in o_indices:
        # Check for ester (acyl)
        ester_matches = mol.GetSubstructMatches(ester_pattern)
        if any(o_idx == match[0] for match in ester_matches):
            substituent_types.append('acyl')
            continue
        
        # Check for vinyl ether (alkenyl)
        vinyl_matches = mol.GetSubstructMatches(vinyl_ether_pattern)
        if any(o_idx == match[0] for match in vinyl_matches):
            substituent_types.append('alkenyl')
            continue
        
        # Default to alkyl (ether)
        substituent_types.append('alkyl')

    # Verify all substituents are valid types
    if all(st in ['acyl', 'alkyl', 'alkenyl'] for st in substituent_types):
        return True, f"Triradylglycerol with substituents: {substituent_types}"
    else:
        return False, f"Invalid substituent types: {substituent_types}"