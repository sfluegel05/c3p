"""
Classifies: CHEBI:24128 furanocoumarin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin (a furan ring fused with a coumarin).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least one 5-membered ring (furan) and one 6-membered ring (coumarin)
    if not any(len(ring) == 5 for ring in rings.AtomRings()):
        return False, "No 5-membered rings (furan) found"
    if not any(len(ring) == 6 for ring in rings.AtomRings()):
        return False, "No 6-membered rings (coumarin) found"

    # Find all furan and coumarin rings
    furan_rings = []
    coumarin_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 5:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if any(atom.GetSymbol() == 'O' for atom in atoms):
                furan_rings.append(ring)
        elif len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if any(atom.GetSymbol() == 'O' for atom in atoms):
                coumarin_rings.append(ring)

    if not furan_rings:
        return False, "No furan rings found"
    if not coumarin_rings:
        return False, "No coumarin rings found"

    # Check if furan and coumarin rings are fused
    for furan_ring in furan_rings:
        for coumarin_ring in coumarin_rings:
            if set(furan_ring).intersection(set(coumarin_ring)):
                return True, "Furanocoumarin structure found"
    
    return False, "No fused furan and coumarin rings found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24128',
                          'name': 'furanocoumarin',
                          'definition': 'Any furochromene that consists of a '
                                        'furan ring fused with a coumarin. The '
                                        'fusion may occur in different ways in '
                                        'give several isomers.',
                          'parents': ['CHEBI:39432']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 20,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}