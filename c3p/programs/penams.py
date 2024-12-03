"""
Classifies: CHEBI:35992 penams
"""
from rdkit import Chem

def is_penams(smiles: str):
    """
    Determines if a molecule is a penam.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penam, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the penam core structure (4-thia-1-azabicyclo[3.2.0]heptan-7-one)
    penam_core_smarts = "C12SC(C)(C)[C@@H](N1C(=O)[C@H]2)C(=O)O"
    penam_core = Chem.MolFromSmarts(penam_core_smarts)
    
    if penam_core is None:
        return False, "Invalid penam core SMARTS"

    # Check if the molecule contains the penam core structure
    if not mol.HasSubstructMatch(penam_core):
        return False, "Molecule does not contain the penam core structure"

    return True, "Molecule contains the penam core structure"

# Example usage
smiles_list = [
    "[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)C(CC)Oc1ccccc1)C(O)=O",  # propicillin
    "[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)[C@]([H])(NC(=O)c1cnc2cccnc2c1O)c1ccccc1)C([O-])=O",  # apalcillin(1-)
    "C=1C=CC=C(C1)CC(N[C@]2([C@@]3(N(C2=O)[C@](C(S3)(C)C)(C(OCOC(=O)C)=O)[H])[H])[H])=O",  # penamecillin
    "[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)Cc1ccccc1)C([O-])=O",  # benzylpenicillin(1-)
    "[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)C(C(=O)Oc1ccc2CCCc2c1)c1ccccc1)C([O-])=O",  # carindacillin(1-)
    "CC1([C@@H](N2[C@@H](S1)[C@@H](C2=O)NC(=O)C(C3=CSC=C3)C(=O)O)C(=O)O)C",  # (2S,5S,6R)-6-[[2-carboxy-1-oxo-2-(3-thiophenyl)ethyl]amino]-3,3-dimethyl-7-oxo-4-thia-1-azabicyclo[3.2.0]heptane-2-carboxylic acid
    "[H]C(=N[C@@H]1C(=O)N2[C@@H](C(=O)OCOC(=O)C(C)(C)C)C(C)(C)S[C@]12[H])N1CCCCCC1",  # pivmecillinam
    "[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)c1c(OCC)ccc2ccccc12)C([O-])=O",  # nafcillin(1-)
    "[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)C(c1ccccc1)S([O-])(=O)=O)C([O-])=O",  # sulbenicillin(2-)
    "[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)C1(N)CCCCC1)C(O)=O",  # cyclacillin
    "[H][C@]1(N[C@@H](C([O-])=O)C(C)(C)S1)[C@H](NC(=O)CC1=CC=CC=C1)C(=O)NCCCC"  # benzylpenicilloyl-butylamine(1-)
]

for smiles in smiles_list:
    result, reason = is_penams(smiles)
    print(f"SMILES: {smiles}")
    print(f"Is penam: {result}")
    print(f"Reason: {reason}")
    print()


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35992',
                          'name': 'penams',
                          'definition': 'Natural and synthetic antibiotics '
                                        'containing the '
                                        '4-thia-1-azabicyclo[3.2.0]heptan-7-one '
                                        'structure, generally assumed to have '
                                        'the 5R configuration unless otherwise '
                                        'specified.',
                          'parents': [   'CHEBI:27171',
                                         'CHEBI:27933',
                                         'CHEBI:38101',
                                         'CHEBI:38106']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: '
              '[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)C(CC)Oc1ccccc1)C(O)=O\n'
              'Is penam: True\n'
              'Reason: Molecule contains the penam core structure\n'
              '\n'
              'SMILES: '
              '[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)[C@]([H])(NC(=O)c1cnc2cccnc2c1O)c1ccccc1)C([O-])=O\n'
              'Is penam: True\n'
              'Reason: Molecule contains the penam core structure\n'
              '\n'
              'SMILES: '
              'C=1C=CC=C(C1)CC(N[C@]2([C@@]3(N(C2=O)[C@](C(S3)(C)C)(C(OCOC(=O)C)=O)[H])[H])[H])=O\n'
              'Is penam: True\n'
              'Reason: Molecule contains the penam core structure\n'
              '\n'
              'SMILES: '
              '[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)Cc1ccccc1)C([O-])=O\n'
              'Is penam: True\n'
              'Reason: Molecule contains the penam core structure\n'
              '\n'
              'SMILES: '
              '[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)C(C(=O)Oc1ccc2CCCc2c1)c1ccccc1)C([O-])=O\n'
              'Is penam: True\n'
              'Reason: Molecule contains the penam core structure\n'
              '\n'
              'SMILES: '
              'CC1([C@@H](N2[C@@H](S1)[C@@H](C2=O)NC(=O)C(C3=CSC=C3)C(=O)O)C(=O)O)C\n'
              'Is penam: True\n'
              'Reason: Molecule contains the penam core structure\n'
              '\n'
              'SMILES: '
              '[H]C(=N[C@@H]1C(=O)N2[C@@H](C(=O)OCOC(=O)C(C)(C)C)C(C)(C)S[C@]12[H])N1CCCCCC1\n'
              'Is penam: True\n'
              'Reason: Molecule contains the penam core structure\n'
              '\n'
              'SMILES: '
              '[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)c1c(OCC)ccc2ccccc12)C([O-])=O\n'
              'Is penam: True\n'
              'Reason: Molecule contains the penam core structure\n'
              '\n'
              'SMILES: '
              '[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)C(c1ccccc1)S([O-])(=O)=O)C([O-])=O\n'
              'Is penam: True\n'
              'Reason: Molecule contains the penam core structure\n'
              '\n'
              'SMILES: '
              '[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)C1(N)CCCCC1)C(O)=O\n'
              'Is penam: True\n'
              'Reason: Molecule contains the penam core structure\n'
              '\n'
              'SMILES: '
              '[H][C@]1(N[C@@H](C([O-])=O)C(C)(C)S1)[C@H](NC(=O)CC1=CC=CC=C1)C(=O)NCCCC\n'
              'Is penam: False\n'
              'Reason: Molecule does not contain the penam core structure\n'
              '\n',
    'num_true_positives': 10,
    'num_false_positives': 0,
    'num_true_negatives': 11,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9090909090909091,
    'f1': 0.9523809523809523,
    'accuracy': None}