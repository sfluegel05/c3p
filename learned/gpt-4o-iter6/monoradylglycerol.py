"""
Classifies: CHEBI:76575 monoradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Look for glycerol backbone pattern (OCC(O)CO)
    glycerol_pattern = Chem.MolFromSmarts("[O][C@@H](CO)[O]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Check for one substituent group - acyl (OC(=O)C), alkyl (OCC), or alk-1-enyl (OC=C)
    acyl_pattern = Chem.MolFromSmarts("OC(=O)[C;!R]")
    alkyl_pattern = Chem.MolFromSmarts("O[C;!R]")
    alkenyl_pattern = Chem.MolFromSmarts("O[C]=C")
    
    acyl_matches = len(mol.GetSubstructMatches(acyl_pattern))
    alkyl_matches = len(mol.GetSubstructMatches(alkyl_pattern))
    alkenyl_matches = len(mol.GetSubstructMatches(alkenyl_pattern))
    
    total_matches = acyl_matches + alkyl_matches + alkenyl_matches
    
    if total_matches != 1:
        return False, f"Expected 1 substituent group, found {total_matches}"

    # Confirm the substituent is a long chain by checking carbon count (typically 14+ for lipids)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 14:
        return False, "Chain too short for typical lipid substituent"
    
    return True, "Contains glycerol backbone with one acyl, alkyl, or alk-1-enyl substituent"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17855',
                          'name': 'monoradylglycerol',
                          'definition': 'Any lipid that is glycerol bearing a single acyl, alkyl or alk-1-enyl substituent at an unspecified position.'},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': ''}