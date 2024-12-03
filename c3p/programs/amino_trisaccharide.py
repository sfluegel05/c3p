"""
Classifies: CHEBI:59266 amino trisaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_amino_trisaccharide(smiles: str):
    """
    Determines if a molecule is an amino trisaccharide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino trisaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of 3 sugar rings
    rings = mol.GetRingInfo()
    sugar_rings = [ring for ring in rings.AtomRings() if len(ring) in [5, 6]]
    if len(sugar_rings) < 3:
        return False, "Less than 3 sugar rings found"

    # Check for the presence of amino groups
    amino_groups = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'N']
    if len(amino_groups) == 0:
        return False, "No amino groups found"

    # Check if amino groups are substituted or unsubstituted
    for atom in amino_groups:
        if atom.GetTotalNumHs() > 0 or any(neighbor.GetSymbol() == 'H' for neighbor in atom.GetNeighbors()):
            return True, "Amino trisaccharide with unsubstituted amino group(s)"
        else:
            return True, "Amino trisaccharide with substituted amino group(s)"

    return None, None  # If the classification cannot be determined

# Example usage:
# smiles = "your_smiles_string_here"
# result, reason = is_amino_trisaccharide(smiles)
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59266',
                          'name': 'amino trisaccharide',
                          'definition': 'An amino oligosaccharide that is a '
                                        'trisaccharide having one or more '
                                        'substituted or unsubstituted amino '
                                        'groups in place of hydroxy groups at '
                                        'unspecified positions.',
                          'parents': ['CHEBI:22483', 'CHEBI:63571']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 28,
    'num_false_positives': 8,
    'num_true_negatives': 12,
    'num_false_negatives': 1,
    'precision': 0.7777777777777778,
    'recall': 0.9655172413793104,
    'f1': 0.8615384615384615,
    'accuracy': None}