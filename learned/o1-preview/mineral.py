"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: CHEBI:46662 mineral
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    A mineral is an inorganic substance, typically formed through geological processes,
    consisting of metals and non-metals, and lacking complex organic structures.

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

    # Check for presence of carbon-carbon bonds
    has_c_c_bond = False
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if atom1.GetSymbol() == 'C' and atom2.GetSymbol() == 'C':
            has_c_c_bond = True
            break
    if has_c_c_bond:
        return False, "Contains carbon-carbon bonds, indicative of organic compounds"

    # Check for presence of rings
    sssr = Chem.GetSSSR(mol)
    if sssr > 0:
        return False, f"Contains {sssr} ring(s), which is not typical for minerals"

    # Check for presence of functional groups typical of organic molecules
    # (e.g., alcohols, amines, alkenes, aromatic rings)
    organic_functional_groups = [
        Chem.MolFromSmarts("[CX3]=[CX3]"),     # Alkenes
        Chem.MolFromSmarts("[CX3]#[CX2]"),     # Alkynes
        Chem.MolFromSmarts("[CX3](=O)[OX2H1]"),# Carboxylic acids
        Chem.MolFromSmarts("[NX3][CX3]=[OX1]"),# Amides
        Chem.MolFromSmarts("[OX2H]"),          # Alcohols
        Chem.MolFromSmarts("[NX3H2]"),         # Primary amines
        Chem.MolFromSmarts("[SX2H]"),          # Thiols
        Chem.MolFromSmarts("c")]               # Aromatic carbons

    for fg in organic_functional_groups:
        if mol.HasSubstructMatch(fg):
            return False, "Contains functional groups typical of organic molecules"

    return True, "Molecule meets criteria for mineral classification"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:46662',
        'name': 'mineral',
        'definition': 'A naturally occurring solid formed through geological processes that has characteristic chemical composition, a highly ordered atomic structure, and specific physical properties.',
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