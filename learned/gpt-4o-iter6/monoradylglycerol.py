"""
Classifies: CHEBI:76575 monoradylglycerol
"""
from rdkit import Chem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol is a glycerol with one acyl, alkyl, or alk-1-enyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Redefine glycerol backbone more flexibly (OCC(O)CO)
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Check for one substituent group - acyl (C=O), alkyl (C), or alk-1-enyl (C=C)
    # More flexible low-key search. Use adjusted smarts for broader match.
    acyl_pattern = Chem.MolFromSmarts("COC(=O)C")
    alkyl_pattern = Chem.MolFromSmarts("COC")
    alkenyl_pattern = Chem.MolFromSmarts("COC=C")
    
    acyl_matches = len(mol.GetSubstructMatches(acyl_pattern))
    alkyl_matches = len(mol.GetSubstructMatches(alkyl_pattern))
    alkenyl_matches = len(mol.GetSubstructMatches(alkenyl_pattern))
    
    # Counting these in an exclusive way could be dependent on found SMILES patterns
    total_matches = acyl_matches + alkyl_matches + alkenyl_matches
    
    if total_matches != 1:
        return False, f"Expected 1 substituent group, found {total_matches} of acyl, alkyl, or alk-1-enyl"

    # Confirm the substituent is in a reasonable chain length for typical lipid (usually > 10 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Chain too short for typical lipid substituent"
    
    return True, "Contains glycerol backbone with one acyl, alkyl, or alk-1-enyl substituent"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17855',
                          'name': 'monoradylglycerol',
                          'definition': 'Any lipid that is glycerol bearing a single acyl, alkyl or alk-1-enyl substituent at an unspecified position.'},
    'message': None,
    'success': True,
    'best': True,
    'error': ''}