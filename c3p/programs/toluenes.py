"""
Classifies: CHEBI:27024 toluenes
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_toluenes(smiles: str):
    """
    Determines if a molecule is a toluene (benzene with one and only one methyl group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a toluene, False otherwise
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

    # Check if all atoms in the aromatic ring are carbon
    for ring in aromatic_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if not all(atom.GetSymbol() == 'C' for atom in atoms):
            return False, "Ring contains non-carbon atoms"

    # Check substituents
    ring_atoms = set(aromatic_rings[0])
    substituents = []

    for atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in ring_atoms:
                substituents.append(neighbor)

    # Check for exactly one methyl group
    methyl_count = 0
    for substituent in substituents:
        if substituent.GetSymbol() == 'C' and substituent.GetTotalNumHs() == 3:
            methyl_count += 1

    if methyl_count == 1:
        return True, "Molecule is a toluene"
    else:
        return False, f"Molecule has {methyl_count} methyl groups instead of one"



__metadata__ = {   'chemical_class': {   'id': 'CHEBI:27024',
                          'name': 'toluenes',
                          'definition': 'Any member of the class of benzenes '
                                        'that is a substituted benzene in '
                                        'which the substituents include one '
                                        '(and only one) methyl group.',
                          'parents': ['CHEBI:22712']},
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
    'num_true_positives': 18,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 5,
    'precision': 0.9473684210526315,
    'recall': 0.782608695652174,
    'f1': 0.8571428571428571,
    'accuracy': None}