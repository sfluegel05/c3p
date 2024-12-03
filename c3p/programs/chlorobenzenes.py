"""
Classifies: CHEBI:23132 chlorobenzenes
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_chlorobenzenes(smiles: str):
    """
    Determines if a molecule is a chlorobenzene (benzene ring substituted by one or more chlorines).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chlorobenzene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered aromatic ring
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

    # Check if at least one aromatic ring is substituted with chlorine
    for ring in aromatic_rings:
        ring_atoms = set(ring)
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'Cl':
                    return True, "Contains chlorobenzene ring"
    
    return False, "No chlorobenzene ring found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23132',
                          'name': 'chlorobenzenes',
                          'definition': 'Any organochlorine compound '
                                        'containing a benzene ring which is '
                                        'substituted by one or more chlorines.',
                          'parents': ['CHEBI:22712', 'CHEBI:36683']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 26-27: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}