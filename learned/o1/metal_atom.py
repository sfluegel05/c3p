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
    # Parse SMILES without sanitization to handle isotopes properly
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return False, "Invalid SMILES string"
    except Exception as e:
        return False, f"Error parsing SMILES: {e}"

    # Attempt to sanitize the molecule; handle valence exceptions
    try:
        Chem.SanitizeMol(mol)
    except Chem.rdchem.AtomValenceException:
        # Ignore valence errors, proceed without full sanitization
        pass
    except Exception as e:
        return False, f"Error during sanitization: {e}"

    # Ensure the molecule contains exactly one atom (allowing for hydrogens)
    if mol.GetNumAtoms() != 1:
        return False, f"Molecule has {mol.GetNumAtoms()} atoms, expected 1"

    # Get the atom
    atom = mol.GetAtomWithIdx(0)
    atomic_num = atom.GetAtomicNum()
    element_symbol = atom.GetSymbol()
    formal_charge = atom.GetFormalCharge()

    # Check if the atom is neutral
    if formal_charge != 0:
        return False, f"Atom is charged: {element_symbol}{formal_charge:+}"

    # Ensure the atom has no bonds (it's not bonded to any other atoms)
    if atom.GetDegree() != 0:
        return False, f"Atom is bonded to other atoms: degree {atom.GetDegree()}"

    # Use RDKit's PeriodicTable to check if the element is a metal
    pt = PeriodicTable.GetPeriodicTable()
    element = pt.GetElementSymbol(atomic_num)

    # List of metallic element categories from the periodic table
    metallic_classes = [
        'Alkali metal',
        'Alkaline earth metal',
        'Lanthanide',
        'Actinide',
        'Transition metal',
        'Post-transition metal',
    ]

    # Get element properties
    group = pt.GetGroup(atomic_num)
    period = pt.GetPeriod(atomic_num)
    block = pt.GetBlock(atomic_num)
    is_metal = pt.IsMetal(atomic_num)

    # Check if element is a metal
    if is_metal:
        return True, f"Atom is a metal: {element_symbol}"
    else:
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
    'attempt': 3,
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