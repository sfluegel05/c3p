"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: CHEBI:17977 diradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is defined as a glycerol backbone (three-carbon chain) with two acyl, alkyl, or alk-1-enyl linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with two or three oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[C@H](O)[CH](O)[CH2O]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for the two ester groups (-O-C(=O)-) or ether/alkyl groups
    ester_pattern = Chem.MolFromSmarts("[O][CX3](=[O])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    ether_alkyl_pattern = Chem.MolFromSmarts("[O][CH2][C]")
    ether_alkyl_matches = mol.GetSubstructMatches(ether_alkyl_pattern)
    
    if len(ester_matches) + len(ether_alkyl_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups and {len(ether_alkyl_matches)} ether/alkyl groups, need exactly 2 linkages"

    return True, "Contains glycerol backbone with two ester and/or ether/alkyl linkages"