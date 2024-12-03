"""
Classifies: CHEBI:24646 hydroquinones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_hydroquinones(smiles: str):
    """
    Determines if a molecule is a hydroquinone (Benzenediols that have the hydroxy substituents in the 1- and 4-positions).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroquinone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered ring
    if not any(len(ring) == 6 for ring in rings.AtomRings()):
        return False, "No 6-membered rings found"

    # Find all aromatic 6-membered rings
    aromatic_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() for atom in atoms):
                aromatic_rings.append(ring)

    if not aromatic_rings:
        return False, "No aromatic 6-membered rings found"

    # Check for hydroxy substituents in the 1- and 4-positions
    for ring in aromatic_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        hydroxy_positions = []
        for idx, atom in enumerate(atoms):
            if atom.GetSymbol() == 'C':
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O' and neighbor.GetTotalDegree() == 1:
                        hydroxy_positions.append(idx)
        
        if 0 in hydroxy_positions and 3 in hydroxy_positions:
            return True, "Hydroquinone with hydroxy substituents in the 1- and 4-positions"

    return False, "No hydroxy substituents in the 1- and 4-positions found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24646',
                          'name': 'hydroquinones',
                          'definition': 'Benzenediols that have the hydroxy '
                                        'substituents in the 1- and '
                                        '4-positions.',
                          'parents': ['CHEBI:33570']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[20:35:18] SMILES Parse Error: syntax error while parsing: '
             'O=C\x01OC(C2=C(O)C=CC(=C2)O)=C/C1=C/C=C/[C@@]3(O[C@@H](C(O)(C)C)CC3)C\n'
             '[20:35:18] SMILES Parse Error: Failed parsing SMILES '
             "'O=C\x01OC(C2=C(O)C=CC(=C2)O)=C/C1=C/C=C/[C@@]3(O[C@@H](C(O)(C)C)CC3)C' "
             'for input: '
             "'O=C\x01OC(C2=C(O)C=CC(=C2)O)=C/C1=C/C=C/[C@@]3(O[C@@H](C(O)(C)C)CC3)C'\n",
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 17,
    'num_false_negatives': 17,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}