"""
Classifies: CHEBI:35504 addition compound
"""
from rdkit import Chem

def is_addition_compound(smiles: str):
    """
    Determines if a molecule is an addition compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an addition compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains multiple distinct components
    components = Chem.GetMolFrags(mol, asMols=True)
    if len(components) < 2:
        return False, "Molecule does not contain multiple distinct components"

    # Check if the molecule contains any donor-acceptor complexes or lattice compounds
    donor_atoms = ['N', 'O', 'S', 'P']  # Common donor atoms
    acceptor_atoms = ['F', 'Cl', 'Br', 'I', 'O', 'N']  # Common acceptor atoms

    has_donor_acceptor = False
    for component in components:
        donor_found = any(atom.GetSymbol() in donor_atoms for atom in component.GetAtoms())
        acceptor_found = any(atom.GetSymbol() in acceptor_atoms for atom in component.GetAtoms())
        if donor_found and acceptor_found:
            has_donor_acceptor = True
            break

    if has_donor_acceptor:
        return True, "Molecule contains donor-acceptor complexes"

    # Check for lattice compounds (presence of metal ions and water molecules)
    metal_ions = ['Li', 'Na', 'K', 'Rb', 'Cs', 'Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Al', 'Ga', 'In', 'Tl', 'Sn', 'Pb', 'Bi']
    water_smiles = Chem.MolFromSmiles("O")
    
    has_metal_ions = any(atom.GetSymbol() in metal_ions for atom in mol.GetAtoms())
    has_water = any(Chem.MolToSmiles(component) == Chem.MolToSmiles(water_smiles) for component in components)

    if has_metal_ions and has_water:
        return True, "Molecule contains metal ions and water molecules (lattice compound)"

    return False, "Molecule does not meet the criteria for an addition compound"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35504',
                          'name': 'addition compound',
                          'definition': 'An addition compound contains two or '
                                        'more simpler compounds that can be '
                                        'packed in a definite ratio into a '
                                        'crystal. The term covers '
                                        'donor-acceptor complexes (adducts) '
                                        'and a variety of lattice compounds.',
                          'parents': ['CHEBI:37577']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[23:31:39] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[23:31:39] WARNING: not removing hydrogen atom without '
             'neighbors\n',
    'stdout': '',
    'num_true_positives': 21,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}