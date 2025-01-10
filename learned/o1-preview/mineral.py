"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: CHEBI:46662 mineral
"""
from rdkit import Chem

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    A mineral is a naturally occurring chemical substance formed through geological processes,
    typically inorganic, and has a characteristic chemical composition.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # List of elements commonly found in minerals
    mineral_elements = {'H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al',
                        'Si', 'P', 'S', 'Cl', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn',
                        'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br',
                        'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',
                        'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Cs', 'Ba', 'La',
                        'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho',
                        'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir',
                        'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Th', 'U'}

    # Get set of elements in the molecule
    atom_symbols = set(atom.GetSymbol() for atom in mol.GetAtoms())
    for element in atom_symbols:
        if element not in mineral_elements:
            return False, f"Element '{element}' not commonly found in minerals"

    # Remove check for carbon-carbon bonds and rings, as minerals can contain them

    # Check for presence of complex organic functional groups
    # Define SMARTS patterns for functional groups typical of organic molecules
    organic_functional_groups = [
        Chem.MolFromSmarts("[#6][#6][#6][#6]"),  # Long carbon chains (4 or more carbons)
        Chem.MolFromSmarts("[#6]=[#6]"),         # Alkenes
        Chem.MolFromSmarts("[#6]#[#6]"),         # Alkynes
        Chem.MolFromSmarts("[#6][OX2H]"),        # Alcohols
        Chem.MolFromSmarts("[#6][NX3]"),         # Amines
        Chem.MolFromSmarts("c1ccccc1"),          # Benzene ring
        Chem.MolFromSmarts("[#6]=O"),            # Carbonyl groups
        Chem.MolFromSmarts("[#6]C(=O)[#6]"),     # Ketones
        Chem.MolFromSmarts("[#6]C(=O)O[#6]"),    # Esters
        Chem.MolFromSmarts("[#6]C(=O)O[H]"),     # Carboxylic acids
        Chem.MolFromSmarts("[#6][#7][#6]"),      # Secondary amines
        Chem.MolFromSmarts("[#6][#16][#6]"),     # Thioethers
    ]

    for fg in organic_functional_groups:
        if mol.HasSubstructMatch(fg):
            return False, "Contains functional groups typical of organic molecules"

    return True, "Molecule meets criteria for mineral classification"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:46662',
        'name': 'mineral',
        'definition': 'In general, a mineral is a chemical substance that is normally crystalline formed and has been formed as a result of geological processes. The term also includes metamict substances (naturally occurring, formerly crystalline substances whose crystallinity has been destroyed by ionising radiation) and can include naturally occurring amorphous substances that have never been crystalline (\'mineraloids\') such as georgite and calciouranoite as well as substances formed by the action of geological processes on bigenic compounds (\'biogenic minerals\').',
        'parents': ['CHEBI:24437']
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
    'attempt': 1,
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