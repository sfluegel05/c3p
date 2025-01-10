"""
Classifies: CHEBI:76579 triradylglycerol
"""
from rdkit import Chem

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.
    A triradylglycerol has a glycerol backbone with three substituent groups (acyl, alkyl, or alk-1-enyl).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycerol backbone: C-C-C with oxygens attached
    glycerol_pattern = Chem.MolFromSmarts("[C]([O])[C]([O])[C]([O])")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for three ester, ether, or similar linkages
    ether_or_ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    matches = mol.GetSubstructMatches(ether_or_ester_pattern)
    
    if len(matches) != 3:
        return False, f"Found {len(matches)} linkage groups, need exactly 3"

    return True, "Contains glycerol backbone with acceptable triradylglycerol linkages"