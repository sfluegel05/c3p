"""
Classifies: CHEBI:22698 benzaldehydes
"""
from rdkit import Chem

def is_benzaldehydes(smiles: str):
    """
    Determines if a molecule is a benzaldehyde (formyl substituted benzene ring and its substituted derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a benzaldehyde, False otherwise
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

    # Check for formyl group (C=O) attached to the benzene ring
    formyl_found = False
    for ring in aromatic_rings:
        ring_atoms = set(ring)
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and neighbor.GetDegree() == 2:
                    # Check for double bonded oxygen
                    for sub_neighbor in neighbor.GetNeighbors():
                        if sub_neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(neighbor.GetIdx(), sub_neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            formyl_found = True
                            break
                if formyl_found:
                    break
            if formyl_found:
                break
        if formyl_found:
            break

    if not formyl_found:
        return False, "No formyl group (C=O) attached to the benzene ring found"

    return True, "Benzaldehyde or substituted benzaldehyde found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22698',
                          'name': 'benzaldehydes',
                          'definition': 'Any arenecarbaldehyde that consists '
                                        'of a formyl substituted benzene ring '
                                        'and its substituted derivatives '
                                        'thereof.',
                          'parents': ['CHEBI:33855']},
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
    'num_true_positives': 12,
    'num_false_positives': 2,
    'num_true_negatives': 10,
    'num_false_negatives': 0,
    'precision': 0.8571428571428571,
    'recall': 1.0,
    'f1': 0.923076923076923,
    'accuracy': None}