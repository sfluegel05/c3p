"""
Classifies: CHEBI:22480 amino disaccharide
"""
from rdkit import Chem

def is_amino_disaccharide(smiles: str):
    """
    Determines if a molecule is an amino disaccharide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino disaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule has exactly two sugar rings
    rings = mol.GetRingInfo().AtomRings()
    sugar_rings = [ring for ring in rings if len(ring) == 5 or len(ring) == 6]
    if len(sugar_rings) != 2:
        return False, "Molecule does not have exactly two sugar rings"

    # Check if there is at least one nitrogen atom in the molecule
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'N']
    if not nitrogen_atoms:
        return False, "No nitrogen atoms found in the molecule"

    # Check if nitrogen atoms are connected to sugar rings
    for nitrogen in nitrogen_atoms:
        for neighbor in nitrogen.GetNeighbors():
            if any(neighbor.GetIdx() in ring for ring in sugar_rings):
                return True, "Molecule is an amino disaccharide"

    return False, "Nitrogen atoms are not connected to the sugar rings"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22480',
                          'name': 'amino disaccharide',
                          'definition': 'A disaccharide derivative that is a '
                                        'disaccharide having one or more '
                                        'substituted or unsubstituted amino '
                                        'groups in place of hydroxy groups at '
                                        'unspecified positions.',
                          'parents': ['CHEBI:22483', 'CHEBI:63353']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 29,
    'num_false_positives': 5,
    'num_true_negatives': 15,
    'num_false_negatives': 1,
    'precision': 0.8529411764705882,
    'recall': 0.9666666666666667,
    'f1': 0.90625,
    'accuracy': None}