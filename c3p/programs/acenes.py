"""
Classifies: CHEBI:51269 acenes
"""
from rdkit import Chem

def is_acenes(smiles: str):
    """
    Determines if a molecule is an acene (Polycyclic aromatic hydrocarbons consisting of fused benzene rings in a rectilinear arrangement and their substitution derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Get all aromatic rings
    aromatic_rings = []
    for ring in rings.AtomRings():
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if all(atom.GetIsAromatic() for atom in atoms):
            aromatic_rings.append(ring)

    if not aromatic_rings:
        return False, "No aromatic rings found"

    # Check if the rings are fused in a rectilinear arrangement
    def is_fused_in_rectilinear(arrangement):
        for i in range(len(arrangement) - 1):
            if len(set(arrangement[i]) & set(arrangement[i + 1])) != 2:
                return False
        return True

    # Check if the aromatic rings form a rectilinear arrangement of fused benzene rings
    fused_rings = []
    for ring in aromatic_rings:
        if len(ring) == 6 and all(mol.GetAtomWithIdx(i).GetSymbol() == 'C' for i in ring):
            fused_rings.append(ring)

    if len(fused_rings) < 2:
        return False, "Less than two fused benzene rings"

    if not is_fused_in_rectilinear(fused_rings):
        return False, "Benzene rings are not fused in a rectilinear arrangement"

    return True, "Molecule is classified as an acene"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51269',
                          'name': 'acenes',
                          'definition': 'Polycyclic aromatic hydrocarbons '
                                        'consisting of fused benzene rings in '
                                        'a rectilinear arrangement and their '
                                        'substitution derivatives.',
                          'parents': ['CHEBI:33836']},
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
    'num_true_positives': 9,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 100,
    'precision': 1.0,
    'recall': 0.08256880733944955,
    'f1': 0.15254237288135594,
    'accuracy': None}