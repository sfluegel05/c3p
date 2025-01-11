"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
"""
Classifies: 17alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    A 17alpha-hydroxy steroid is a steroid with an alpha-oriented hydroxyl group at the C17 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure molecule has stereochemistry information
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    # Steroid backbone SMARTS with labeled atoms (C17 is labeled as atom 17)
    steroid_smarts = """
    [#6]-1=[#6]-[#6]-[#6]-2-[#6]-[#6]-3=[#6]-[#6]-[#6]-[#6]-[#6]-3-
    [#6]-[#6]-2-[#6]-1
    """

    # Create a template steroid molecule with labeled atoms
    steroid_template = Chem.MolFromSmiles('C1CCC2C3CCC4=CC(=O)CC[C@]4(C)[C@@H]3CC[C@]12C')
    if steroid_template is None:
        return False, "Failed to create steroid template"

    # Label the C17 atom in the template (atom index 16 in zero-based indexing)
    atom_map = {atom.GetIdx(): atom.GetIdx() + 1 for atom in steroid_template.GetAtoms()}
    steroid_template.GetAtomWithIdx(16).SetAtomMapNum(17)  # C17

    # Perform substructure match to align the molecule to the steroid template
    match = mol.GetSubstructMatch(steroid_template)
    if not match:
        return False, "Molecule does not match steroid backbone"

    # Create an atom map from template to molecule
    atom_map = {template_idx: mol_idx for template_idx, mol_idx in enumerate(match)}

    # Get the C17 atom in the molecule
    c17_idx = atom_map.get(16)  # Template atom index 16 corresponds to C17
    if c17_idx is None:
        return False, "C17 atom not found in molecule"

    c17_atom = mol.GetAtomWithIdx(c17_idx)

    # Check for hydroxyl group at C17
    has_oh = False
    for neighbor in c17_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 8:
            # Check if the oxygen is part of a hydroxyl group
            if len(neighbor.GetNeighbors()) == 1:
                has_oh = True
                oxygen_idx = neighbor.GetIdx()
                break

    if not has_oh:
        return False, "No hydroxyl group attached to C17"

    # Check stereochemistry at C17 (alpha orientation)
    stereo = c17_atom.GetChiralTag()
    if stereo == Chem.CHI_UNSPECIFIED:
        return False, "Chirality at C17 is unspecified"

    # Determine if the configuration is alpha at C17
    # For steroids, alpha orientation corresponds to S configuration at C17
    # RDKit uses CIP rules to assign R/S configuration
    c17_symbol = c17_atom.GetSymbol()
    if c17_atom.HasProp('_CIPCode'):
        cip_code = c17_atom.GetProp('_CIPCode')
        if cip_code != 'S':
            return False, f"Hydroxyl group at C17 is not alpha-oriented (CIP code: {cip_code})"
    else:
        return False, "Unable to assign CIP code to C17 atom"

    return True, "Molecule is a 17alpha-hydroxy steroid with correct stereochemistry at C17"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35342',
                          'name': '17alpha-hydroxy steroid',
                          'definition': 'The alpha-stereoisomer of 17-hydroxy '
                                        'steroid.',
                          'parents': ['CHEBI:36838'],
                          'xrefs': ['KEGG:C03336'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               '[H][C@@]12CC[C@](O)(C(=O)CO)[C@@]1(C)CC(=O)[C@@]1([H])[C@@]2([H])C[C@H](O)C2=CC(=O)C=C[C@]12C '
               'NAME: 6alpha-hydroxyprednisone REASON: MISSED No tertiary '
               'chiral carbon with alpha-hydroxyl group found\n'
               ' * SMILES: '
               '[H][C@@]12CC[C@@]3([H])[C@]4([H])CC[C@](O)(C(=O)CO)[C@@]4(C)CC(=O)[C@]3([H])[C@@]1(C)CCC(=O)C2 '
               'NAME: 4,5alpha-dihydrocortisone REASON: MISSED No chiral '
               'carbon with attached hydroxyl group found\n'
               ' * SMILES: '
               '[H][C@@]12CCC3CC(=O)CC[C@]3(C)[C@@]1([H])C(=O)C[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(=O)CO '
               'NAME: 4,5-dihydrocortisone REASON: MISSED No chiral carbon '
               'with attached hydroxyl group found\n'
               ' * SMILES: '
               '[H][C@@]12C[C@@H](O)[C@](O)(C(=O)CO)[C@@]1(C)C[C@H](O)[C@@]1([H])[C@@]2([H])CCC2=CC(=O)C=C[C@]12C '
               'NAME: 16alpha-hydroxyprednisolone REASON: MISSED No tertiary '
               'chiral carbon with alpha-hydroxyl group found\n'
               ' * SMILES: '
               '[H][C@@]12C[C@@H](C)[C@](O)(C(=O)CCl)[C@@]1(C)C[C@H](O)[C@@]1(Cl)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C '
               'NAME: mometasone REASON: MISSED No tertiary chiral carbon with '
               'alpha-hydroxyl group found\n'
               ' * SMILES: '
               '[H][C@@]12C[C@H](C)[C@](O)(C(=O)COP(O)(O)=O)[C@@]1(C)C[C@H](O)[C@@]1(F)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C '
               'NAME: betamethasone phosphate REASON: MISSED No tertiary '
               'chiral carbon with alpha-hydroxyl group found\n'
               ' * SMILES: '
               'C1[C@@]2([C@@]([C@]3(CCC(C=C3[C@@H]1O)=O)C)(CC[C@]4([C@]2(C[C@H]([C@]4(O)C(C)=O)O)[H])C)[H])[H] '
               'NAME: 6beta,16alpha,17alpha-trihydroxyprogesterone REASON: '
               'MISSED No tertiary chiral carbon with alpha-hydroxyl group '
               'found\n'
               ' * SMILES: '
               '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@@]1(C)[C@@]2([H])CC[C@]1(O)[C@H](C)O '
               'NAME: (20S)-17,20-dihydroxypregn-4-en-3-one REASON: MISSED No '
               'tertiary chiral carbon with alpha-hydroxyl group found\n'
               ' * SMILES: '
               '[H][C@@]12C[C@@H](O)[C@](O)(C(=O)CO)[C@@]1(C)C[C@H](O)[C@@]1(F)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C '
               'NAME: triamcinolone REASON: MISSED No tertiary chiral carbon '
               'with alpha-hydroxyl group found\n'
               ' * SMILES: '
               '[H][C@@]12CC[C@](O)(C(=O)COC(=O)CCC(O)=O)[C@@]1(C)C[C@H](O)[C@@]1([H])[C@@]2([H])CCC2=CC(=O)C=C[C@]12C '
               'NAME: prednisolone succinate REASON: MISSED No tertiary chiral '
               'carbon with alpha-hydroxyl group found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'CNC(=O)CC[C@H](N)C(O)=O',
                                     'name': 'N(5)-methyl-L-glutamine',
                                     'reason': 'Molecule does not match '
                                               'steroid backbone'},
                                 {   'smiles': 'C[C@H](OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c1nc(N)[nH]c2=O)C(O)=O',
                                     'name': "L-lactyl-2-diphospho-5'-guanosine",
                                     'reason': 'Molecule does not match '
                                               'steroid backbone'},
                                 {   'smiles': '[H][C@@]1(CC[C@]2(C)[C@]1([H])CC[C@]1([H])[C@@]3(C)CCCC(C)(C)[C@]3([H])CC[C@@]21C)C(=C)CCC=C(C)C',
                                     'name': 'dammara-20,24-diene',
                                     'reason': 'Molecule does not match '
                                               'steroid backbone'},
                                 {   'smiles': 'O=C1C2=C(O)C3=C(O)C=CC=C3C(=C2C(C(=O)C)C(C1)(O)C)C=4C=5C(C(O)=C6C4C(C(=O)C)C(O)(C)CC6=O)=C(O)C=CC5',
                                     'name': 'A 39183A',
                                     'reason': 'Molecule does not match '
                                               'steroid backbone'},
                                 {   'smiles': 'C[C@H](O)[C@@H](O)C1=Nc2c(NC1)[nH]c(N)nc2=O',
                                     'name': 'L-threo-7,8-dihydrobiopterin',
                                     'reason': 'Molecule does not match '
                                               'steroid backbone'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@H]1OC)C)C)CC3=CC=C(C=C3)F',
                                     'name': 'N-[(4S,7R,8S)-5-[(4-fluorophenyl)methyl]-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'Molecule does not match '
                                               'steroid backbone'},
                                 {   'smiles': 'C(CCN1CCCCC1)(CCN(C(C)C)C(C)=O)(C(N)=O)C2=C(Cl)C=CC=C2',
                                     'name': 'bidisomide',
                                     'reason': 'Molecule does not match '
                                               'steroid backbone'},
                                 {   'smiles': 'OCCCCCCCCCCC#CC=C',
                                     'name': '13-Tetradece-11-yn-1-ol',
                                     'reason': 'Molecule does not match '
                                               'steroid backbone'},
                                 {   'smiles': 'C[C@@H]1CCCCO[C@@H]([C@@H](CN(C(=O)C2=C(O1)C=CC(=C2)N(C)C)[C@@H](C)CO)C)CN(C)CC3=CC=C(C=C3)OC',
                                     'name': '(3R,9S,10R)-16-(dimethylamino)-12-[(2S)-1-hydroxypropan-2-yl]-9-[[(4-methoxyphenyl)methyl-methylamino]methyl]-3,10-dimethyl-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-13-one',
                                     'reason': 'Molecule does not match '
                                               'steroid backbone'},
                                 {   'smiles': 'O(CCC12CC3CC(C1)CC(C2)C3)C(=O)CC4=CC=C(OCC(O)CNC(C)C)C=C4',
                                     'name': 'Adaprolol',
                                     'reason': 'Molecule does not match '
                                               'steroid backbone'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@@]12CC[C@](O)(C(=O)CO)[C@@]1(C)CC(=O)[C@@]1([H])[C@@]2([H])C[C@H](O)C2=CC(=O)C=C[C@]12C',
                                      'name': '6alpha-hydroxyprednisone',
                                      'reason': 'Molecule does not match '
                                                'steroid backbone'},
                                  {   'smiles': '[H][C@@]12CC[C@@]3([H])[C@]4([H])CC[C@](O)(C(=O)CO)[C@@]4(C)CC(=O)[C@]3([H])[C@@]1(C)CCC(=O)C2',
                                      'name': '4,5alpha-dihydrocortisone',
                                      'reason': 'Molecule does not match '
                                                'steroid backbone'},
                                  {   'smiles': '[H][C@@]12CCC3CC(=O)CC[C@]3(C)[C@@]1([H])C(=O)C[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(=O)CO',
                                      'name': '4,5-dihydrocortisone',
                                      'reason': 'Molecule does not match '
                                                'steroid backbone'},
                                  {   'smiles': '[H][C@@]12C[C@@H](O)[C@](O)(C(=O)CO)[C@@]1(C)C[C@H](O)[C@@]1([H])[C@@]2([H])CCC2=CC(=O)C=C[C@]12C',
                                      'name': '16alpha-hydroxyprednisolone',
                                      'reason': 'Molecule does not match '
                                                'steroid backbone'},
                                  {   'smiles': '[H][C@@]12C[C@@H](C)[C@](O)(C(=O)CCl)[C@@]1(C)C[C@H](O)[C@@]1(Cl)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C',
                                      'name': 'mometasone',
                                      'reason': 'Molecule does not match '
                                                'steroid backbone'},
                                  {   'smiles': '[H][C@@]12C[C@H](C)[C@](O)(C(=O)COP(O)(O)=O)[C@@]1(C)C[C@H](O)[C@@]1(F)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C',
                                      'name': 'betamethasone phosphate',
                                      'reason': 'Molecule does not match '
                                                'steroid backbone'},
                                  {   'smiles': 'C1[C@@]2([C@@]([C@]3(CCC(C=C3[C@@H]1O)=O)C)(CC[C@]4([C@]2(C[C@H]([C@]4(O)C(C)=O)O)[H])C)[H])[H]',
                                      'name': '6beta,16alpha,17alpha-trihydroxyprogesterone',
                                      'reason': 'No hydroxyl group attached to '
                                                'C17'},
                                  {   'smiles': '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@@]1(C)[C@@]2([H])CC[C@]1(O)[C@H](C)O',
                                      'name': '(20S)-17,20-dihydroxypregn-4-en-3-one',
                                      'reason': 'No hydroxyl group attached to '
                                                'C17'},
                                  {   'smiles': '[H][C@@]12C[C@@H](O)[C@](O)(C(=O)CO)[C@@]1(C)C[C@H](O)[C@@]1(F)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C',
                                      'name': 'triamcinolone',
                                      'reason': 'Molecule does not match '
                                                'steroid backbone'},
                                  {   'smiles': '[H][C@@]12CC[C@](O)(C(=O)COC(=O)CCC(O)=O)[C@@]1(C)C[C@H](O)[C@@]1([H])[C@@]2([H])CCC2=CC(=O)C=C[C@]12C',
                                      'name': 'prednisolone succinate',
                                      'reason': 'Molecule does not match '
                                                'steroid backbone'}],
    'attempt': 1,
    'reasoning': 'The previous code did not correctly identify the '
                 '17alpha-hydroxy steroids because it relied on a general '
                 'SMARTS pattern to find a chiral tertiary carbon with a '
                 'hydroxyl group (`[C@@H](O)[C]`). This pattern may not '
                 'accurately capture the specific structure and '
                 'stereochemistry of the 17alpha-hydroxyl group in steroids. '
                 'Additionally, the code did not map the atoms to the standard '
                 'steroid numbering system, making it difficult to '
                 'specifically identify the C17 position.\n'
                 '\n'
                 'To improve the code, we need to:\n'
                 '\n'
                 '1. **Map the steroid backbone with assigned atom numbers**: '
                 'Use a steroid template with numbered atoms to accurately '
                 'identify the C17 position.\n'
                 '2. **Check for a hydroxyl group at the C17 position**: '
                 'Ensure that there is an -OH group attached to the C17 '
                 'carbon.\n'
                 '3. **Verify the stereochemistry at C17**: Confirm that the '
                 'hydroxyl group is in the alpha orientation at C17.\n'
                 '4. **Ensure correct stereochemistry is assigned**: Use RDKit '
                 'functions to assign stereochemistry to the molecule.\n'
                 '\n'
                 'By implementing these steps, we can more accurately classify '
                 '17alpha-hydroxy steroids.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 13,
    'num_false_positives': 50,
    'num_true_negatives': 142164,
    'num_false_negatives': 73,
    'num_negatives': None,
    'precision': 0.20634920634920634,
    'recall': 0.1511627906976744,
    'f1': 0.174496644295302,
    'accuracy': 0.9991356289529164,
    'negative_predictive_value': 0.9994867720775888}