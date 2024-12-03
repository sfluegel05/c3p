"""
Classifies: CHEBI:33240 coordination entity
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_coordination_entity(smiles: str):
    """
    Determines if a molecule is a coordination entity (an assembly consisting of a central atom, usually metallic, to which is attached a surrounding array of other groups of atoms (ligands)).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a coordination entity, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of a metal atom
    metals = ["Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]
    metal_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() in metals]

    if not metal_atoms:
        return False, "No metal atoms found"

    # Check for ligands attached to metal atom
    for metal in metal_atoms:
        neighbors = metal.GetNeighbors()
        if len(neighbors) > 0:
            return True, "Coordination entity found with metal atom and ligands"
    
    return False, "Metal atom found but no ligands attached"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33240',
                          'name': 'coordination entity',
                          'definition': 'An assembly consisting of a central '
                                        'atom (usually metallic) to which is '
                                        'attached a surrounding array of other '
                                        'groups of atoms (ligands).',
                          'parents': ['CHEBI:37577']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 36-37: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}