"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: CHEBI:60324 hemiaminal
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule contains a hemiaminal group based on its SMILES string.
    A hemiaminal has both an amino group and a hydroxy group attached to the same carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a hemiaminal group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to handle implicit H cases
    mol = Chem.AddHs(mol)

    # Define SMARTS patterns for hemiaminal groups
    # Pattern 1: Carbon with both OH and NH/NH2 attached
    # [$(C([OH1])[NH1,NH2])] matches carbon with both OH and NH/NH2
    hemiaminal_pattern1 = Chem.MolFromSmarts('[$(C([OH1])[NH1,NH2])]')
    
    # Pattern 2: Carbon with both OH and N in ring (covers cyclic cases)
    # [$(C([OH1])([#7D2,#7D3]))] matches carbon with OH and ring nitrogen
    hemiaminal_pattern2 = Chem.MolFromSmarts('[$(C([OH1])([#7D2,#7D3]))]')
    
    # Check for matches
    matches1 = mol.GetSubstructMatches(hemiaminal_pattern1)
    matches2 = mol.GetSubstructMatches(hemiaminal_pattern2)
    
    total_matches = len(matches1) + len(matches2)
    
    if total_matches == 0:
        return False, "No hemiaminal group found (no carbon with both OH and NH/NH2 groups)"

    # Additional validation to avoid false positives
    # Check if any match is part of a carboxylic acid or similar groups
    false_positive_pattern = Chem.MolFromSmarts('[$(C(=O)[OH1])]')  # carboxylic acid pattern
    false_positives = mol.GetSubstructMatches(false_positive_pattern)
    
    # Convert matches to sets for comparison
    matches_set = set([m[0] for m in matches1] + [m[0] for m in matches2])
    false_pos_set = set([m[0] for m in false_positives])
    
    # Remove false positives
    true_matches = matches_set - false_pos_set
    
    if not true_matches:
        return False, "Found potential matches but they are part of carboxylic acid groups"

    reason = f"Found {len(true_matches)} hemiaminal group(s) (carbon atom(s) with both OH and N-containing groups attached)"
    return True, reason