"""
Classifies: CHEBI:83403 monochlorobenzenes
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_monochlorobenzenes(smiles: str):
    """
    Determines if a molecule is a monochlorobenzene (benzene ring with exactly one chlorine substituent).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monochlorobenzene, False otherwise
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

    # Check if any of the aromatic rings have exactly one chlorine substituent
    for ring in aromatic_rings:
        ring_atoms = set(ring)
        chlorine_count = 0
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'Cl' and neighbor.GetIdx() not in ring_atoms:
                    chlorine_count += 1

        if chlorine_count == 1:
            return True, "Monochlorobenzene found"

    return False, "No monochlorobenzene found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:83403',
                          'name': 'monochlorobenzenes',
                          'definition': 'Any member of the class of '
                                        'chlorobenzenes containing a mono- or '
                                        'poly-substituted benzene ring in '
                                        'which only one substituent is '
                                        'chlorine.',
                          'parents': ['CHEBI:23132']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 15-16: malformed \\N character escape (<string>, line '
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