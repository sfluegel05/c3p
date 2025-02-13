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
                        'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'Th', 'U'}

    # Get set of elements in the molecule
    atom_symbols = set(atom.GetSymbol() for atom in mol.GetAtoms())
    for element in atom_symbols:
        if element not in mineral_elements:
            return False, f"Element '{element}' not commonly found in minerals"

    # Exclude molecules with direct metal-carbon bonds (organometallics)
    metals_atomic_nums = [3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24,
                          25, 26, 27, 28, 29, 30, 31, 37, 38, 39,
                          40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
                          55, 56, 57, 58, 59, 60, 61, 62, 63, 64,
                          65, 66, 67, 68, 69, 70, 71, 72, 73, 74,
                          75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
                          90, 92]  # Atomic numbers of metals
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        if ((begin_atom.GetAtomicNum() in metals_atomic_nums and end_atom.GetAtomicNum() == 6) or
            (end_atom.GetAtomicNum() in metals_atomic_nums and begin_atom.GetAtomicNum() == 6)):
            return False, "Contains metal-carbon bond typical of organometallic compounds"

    # Exclude molecules with aromatic rings
    aromatic_ring = Chem.MolFromSmarts("a1aaaaa1")
    if mol.HasSubstructMatch(aromatic_ring):
        return False, "Contains aromatic ring typical of organic molecules"

    # Exclude molecules with long carbon chains (more than 4 contiguous carbons)
    long_carbon_chain = Chem.MolFromSmarts("[CH2]([CH2])[CH2][CH2]")
    if mol.HasSubstructMatch(long_carbon_chain):
        return False, "Contains long carbon chain typical of organic molecules"

    # Exclude molecules with functional groups typical of complex organic molecules
    organic_functional_groups = [
        Chem.MolFromSmarts("[#6]=[#6]"),              # Alkenes
        Chem.MolFromSmarts("[#6]#[#6]"),              # Alkynes
        Chem.MolFromSmarts("[#6][OX2H]"),             # Alcohols
        Chem.MolFromSmarts("[#6][NX3;H2,H1;!$([N][O])]"),  # Amines (primary and secondary)
        Chem.MolFromSmarts("[#6]C(=O)[#6]"),          # Ketones
        Chem.MolFromSmarts("[#6]C(=O)O[#6]"),         # Esters
        Chem.MolFromSmarts("[#6][#16][#6]"),          # Thioethers
        Chem.MolFromSmarts("c"),                      # Any aromatic carbon atom
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
    'attempt': 2,
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