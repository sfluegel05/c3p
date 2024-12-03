"""
Classifies: CHEBI:37667 sesquiterpene lactone
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_sesquiterpene_lactone(smiles: str):
    """
    Determines if a molecule is a sesquiterpene lactone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpene lactone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a lactone ring (cyclic ester)
    lactone_found = False
    for ring in mol.GetRingInfo().BondRings():
        if len(ring) == 5 or len(ring) == 6:
            ring_atoms = [mol.GetBondWithIdx(bond_idx).GetBeginAtomIdx() for bond_idx in ring]
            ring_atoms += [mol.GetBondWithIdx(bond_idx).GetEndAtomIdx() for bond_idx in ring]
            ring_atoms = set(ring_atoms)
            if any(mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'O' for atom_idx in ring_atoms):
                lactone_found = True
                break

    if not lactone_found:
        return False, "No lactone ring found"

    # Check for sesquiterpene structure (15 carbon atoms)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if num_carbons < 15:
        return False, "Not enough carbon atoms for a sesquiterpene structure"

    # Check for the presence of multiple rings (multicyclic structure)
    if len(mol.GetRingInfo().AtomRings()) < 2:
        return False, "Not enough rings for a multicyclic structure"

    return True, "Molecule is a sesquiterpene lactone"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37667',
                          'name': 'sesquiterpene lactone',
                          'definition': 'Any member of a diverse class of '
                                        'complex, multicyclic phytochemicals '
                                        'showing a variety of skeleton '
                                        'arrangements and bioactivities, and '
                                        'having in common a sesquiterpenoid '
                                        'structure including a lactone ring.',
                          'parents': ['CHEBI:26658', 'CHEBI:37668']},
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
    'num_true_positives': 36,
    'num_false_positives': 3,
    'num_true_negatives': 17,
    'num_false_negatives': 0,
    'precision': 0.9230769230769231,
    'recall': 1.0,
    'f1': 0.9600000000000001,
    'accuracy': None}