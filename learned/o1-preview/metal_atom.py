"""
Classifies: CHEBI:33521 metal atom
"""
"""
Classifies: metal atom
"""
from rdkit import Chem
from rdkit.Chem import PeriodicTable

def is_metal_atom(smiles: str):
    """
    Determines if a molecule is a metal atom based on its SMILES string.
    A metal atom is an atom of an element that exhibits typical metallic properties,
    being typically shiny, with high electrical and thermal conductivity.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a metal atom, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule has exactly one atom
    if mol.GetNumAtoms() != 1:
        return False, f"Molecule has {mol.GetNumAtoms()} atoms, expected 1"

    # Get the atom
    atom = mol.GetAtomWithIdx(0)
    atomic_num = atom.GetAtomicNum()

    # List of atomic numbers of metal elements
    metal_atomic_numbers = set([
        # Alkali metals
        3,   # Lithium (Li)
        11,  # Sodium (Na)
        19,  # Potassium (K)
        37,  # Rubidium (Rb)
        55,  # Cesium (Cs)
        87,  # Francium (Fr)
        # Alkaline earth metals
        4,   # Beryllium (Be)
        12,  # Magnesium (Mg)
        20,  # Calcium (Ca)
        38,  # Strontium (Sr)
        56,  # Barium (Ba)
        88,  # Radium (Ra)
        # Lanthanides
        57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
        # Actinides
        89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103,
        # Transition metals
        21, 22, 23, 24, 25, 26, 27, 28, 29, 30,  # Sc to Zn
        39, 40, 41, 42, 43, 44, 45, 46, 47, 48,  # Y to Cd
        72, 73, 74, 75, 76, 77, 78, 79, 80,      # Hf to Hg
        104,105,106,107,108,109,110,111,112,     # Rf to Cn
        # Post-transition metals
        13, 31, 49, 50, 81, 82, 83, 113, 114, 115, 116, # Al, Ga, In, Sn, Tl, Pb, Bi, Nh, Fl, Mc, Lv
        # Metalloids (some may not consider these metals)
        5, 14, 32, 33, 51, 52, 84, # B, Si, Ge, As, Sb, Te, Po
        # Others
        56, 87, 88, # Ba, Fr, Ra
    ])

    # Check if the atomic number is in the set of metals
    if atomic_num in metal_atomic_numbers:
        element_symbol = atom.GetSymbol()
        return True, f"Atom is a metal: {element_symbol}"
    else:
        element_symbol = atom.GetSymbol()
        return False, f"Atom is not a metal: {element_symbol}"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'metal atom',
        'definition': 'An atom of an element that exhibits typical metallic properties, being typically shiny, with high electrical and thermal conductivity.',
        'parents': []
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}