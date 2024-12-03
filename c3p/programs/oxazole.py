"""
Classifies: CHEBI:35790 oxazole
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_oxazole(smiles: str):
    """
    Determines if a molecule is an oxazole (a five-membered heterocyclic aromatic skeleton containing one N and one O atom).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxazole, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least one 5-membered ring
    if not any(len(ring) == 5 for ring in rings.AtomRings()):
        return False, "No 5-membered rings found"

    # Find all 5-membered rings
    five_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 5]

    for ring in five_membered_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        atom_symbols = [atom.GetSymbol() for atom in atoms]
        if 'N' in atom_symbols and 'O' in atom_symbols:
            n_count = sum(1 for atom in atoms if atom.GetSymbol() == 'N')
            o_count = sum(1 for atom in atoms if atom.GetSymbol() == 'O')
            if n_count == 1 and o_count == 1 and all(atom.GetIsAromatic() for atom in atoms):
                return True, "Oxazole ring found"

    return False, "No oxazole ring found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35790',
                          'name': 'oxazole',
                          'definition': 'An azole based on a five-membered '
                                        'heterocyclic aromatic skeleton '
                                        'containing one N and one O atom.',
                          'parents': ['CHEBI:38104', 'CHEBI:68452']},
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
    'num_true_positives': 25,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 5,
    'precision': 1.0,
    'recall': 0.8333333333333334,
    'f1': 0.9090909090909091,
    'accuracy': None}