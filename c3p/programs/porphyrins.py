"""
Classifies: CHEBI:26214 porphyrins
"""
from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for exactly four 5-membered rings (pyrrole rings)
    pyrrole_rings = [ring for ring in rings.AtomRings() if len(ring) == 5]
    if len(pyrrole_rings) != 4:
        return False, "Does not contain exactly four 5-membered rings"

    # Check if each 5-membered ring is a pyrrole ring (contains one nitrogen)
    for ring in pyrrole_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if sum(1 for atom in atoms if atom.GetSymbol() == 'N') != 1:
            return False, "One or more 5-membered rings is not a pyrrole ring"

    # Check if the pyrrole rings are connected through the alpha-positions by methine groups
    # Methine group: -CH-
    alpha_positions = []
    for ring in pyrrole_rings:
        # Find the nitrogen atom in the ring
        nitrogen_idx = next(i for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == 'N')
        # Alpha positions are the carbon atoms adjacent to the nitrogen
        alpha_positions.append([i for i in ring if mol.GetBondBetweenAtoms(nitrogen_idx, i) is not None and mol.GetAtomWithIdx(i).GetSymbol() == 'C'])

    # Check if the alpha positions are connected by methine groups
    connected_by_methine = True
    for i in range(len(alpha_positions)):
        current_alpha = alpha_positions[i]
        next_alpha = alpha_positions[(i + 1) % len(alpha_positions)]
        if not any(mol.GetBondBetweenAtoms(a, b) is not None and mol.GetAtomWithIdx(b).GetSymbol() == 'C' and mol.GetAtomWithIdx(b).GetDegree() == 3 for a in current_alpha for b in next_alpha):
            connected_by_methine = False
            break

    if not connected_by_methine:
        return False, "Pyrrole rings are not connected through alpha-positions by methine groups"

    return True, "Molecule is a porphyrin"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26214',
                          'name': 'porphyrins',
                          'definition': 'Natural pigments containing a '
                                        'fundamental skeleton of four pyrrole '
                                        'nuclei united through the '
                                        'alpha-positions by four methine '
                                        'groups to form a macrocyclic '
                                        'structure.',
                          'parents': ['CHEBI:36309']},
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
    'num_true_negatives': 16,
    'num_false_negatives': 16,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}