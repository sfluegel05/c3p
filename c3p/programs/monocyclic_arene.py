"""
Classifies: CHEBI:33847 monocyclic arene
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_monocyclic_arene(smiles: str):
    """
    Determines if a molecule is a monocyclic arene (monocyclic aromatic hydrocarbon).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monocyclic arene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for exactly one aromatic ring
    aromatic_rings = []
    for ring in rings.AtomRings():
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if all(atom.GetIsAromatic() for atom in atoms):
            aromatic_rings.append(ring)

    if len(aromatic_rings) != 1:
        return False, "Molecule does not have exactly one aromatic ring"

    # Check if all atoms in the aromatic ring are carbon
    ring = aromatic_rings[0]
    atoms = [mol.GetAtomWithIdx(i) for i in ring]
    if not all(atom.GetSymbol() == 'C' for atom in atoms):
        return False, "Aromatic ring contains non-carbon atoms"

    return True, "Molecule is a monocyclic arene"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33847',
                          'name': 'monocyclic arene',
                          'definition': 'A monocyclic aromatic hydrocarbon.',
                          'parents': ['CHEBI:33658']},
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
    'num_true_positives': 21,
    'num_false_positives': 4,
    'num_true_negatives': 16,
    'num_false_negatives': 8,
    'precision': 0.84,
    'recall': 0.7241379310344828,
    'f1': 0.7777777777777777,
    'accuracy': None}