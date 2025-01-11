"""
Classifies: CHEBI:50128 biflavonoid
"""
"""
Classifies: biflavonoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    A biflavonoid is a flavonoid dimer obtained by the coupling of at least two flavonoid units (aryl-substituted benzopyran rings or derivatives), resulting in two ring systems being joined together by a single atom or bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a biflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general flavonoid core SMARTS pattern (C6-C3-C6 skeleton)
    flavonoid_core_smarts = '[#6]1~[#6]~[#6]~[#6]~[#6]~[#6]1~[#6]2~[#6]~[#6]~[#6]~[#6]~[#8]~2'  # General flavonoid core
    flavonoid_core = Chem.MolFromSmarts(flavonoid_core_smarts)
    if flavonoid_core is None:
        return False, "Invalid flavonoid core SMARTS pattern"

    # Find matches for the flavonoid core
    flavonoid_matches = mol.GetSubstructMatches(flavonoid_core)
    num_flavonoid_cores = len(flavonoid_matches)

    if num_flavonoid_cores < 2:
        return False, f"Found {num_flavonoid_cores} flavonoid core(s), need at least 2"

    # Create sets of atoms for each flavonoid core
    flavonoid_core_atoms = [set(match) for match in flavonoid_matches]

    # Check if flavonoid cores are connected via a single atom or bond
    found_connected = False
    for i in range(num_flavonoid_cores):
        for j in range(i+1, num_flavonoid_cores):
            # Check for direct bonds or shared atoms between cores i and j
            inter_core_bonds = 0
            inter_core_atoms = flavonoid_core_atoms[i] & flavonoid_core_atoms[j]
            if inter_core_atoms:
                inter_core_bonds = len(inter_core_atoms)
            else:
                for atom_i in flavonoid_core_atoms[i]:
                    for atom_j in flavonoid_core_atoms[j]:
                        bond = mol.GetBondBetweenAtoms(atom_i, atom_j)
                        if bond is not None:
                            inter_core_bonds += 1
            if inter_core_bonds == 1:
                found_connected = True
                break
        if found_connected:
            break

    if not found_connected:
        return False, "Flavonoid cores are not connected via a single atom or bond"

    return True, "Contains at least two flavonoid cores connected via a single atom or bond"


__metadata__ = {   'chemical_class': {   'name': 'biflavonoid',
                              'definition': 'A flavonoid oligomer that is obtained by the oxidative coupling of at least two units of aryl-substituted benzopyran rings or its substituted derivatives, resulting in the two ring systems being joined together by a single atom or bond.'},
        'config': {   'llm_model_name': 'lbl/claude-sonnet'},
        'success': True}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50128',
                          'name': 'biflavonoid',
                          'definition': 'A flavonoid oligomer that is obtained '
                                        'by the oxidative coupling of at least '
                                        'two units of aryl-substituted '
                                        'benzopyran rings or its substituted '
                                        'derivatives, resulting in the two '
                                        'ring systems being joined together by '
                                        'a single atom or bond.',
                          'parents': ['CHEBI:72720'],
                          'xrefs': ['Wikipedia:Biflavonoid'],
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
               'COc1cc2OC(CC(=O)c2c(O)c1Oc1ccc(cc1)-c1cc(=O)c2c(O)c(C)c(O)cc2o1)c1ccc(O)cc1 '
               'NAME: taiwanhomoflavone-B REASON: MISSED Found 0 flavonoid '
               'core(s), need at least 2\n'
               ' * SMILES: '
               '[H][C@]1(Oc2cc(O)cc(O)c2C(=O)[C@@]1([H])[C@@]1([H])C(=O)c2c(O)cc(O)cc2O[C@]1([H])c1ccc(OC)cc1)c1ccc(OC)cc1 '
               'NAME: chamaejasmenin A REASON: MISSED Found 0 flavonoid '
               'core(s), need at least 2\n'
               ' * SMILES: '
               'O[C@H]1Cc2c(O)cc(O)c([C@H]3[C@H](O)[C@H](Oc4cc(O)cc(O)c34)c3ccc(O)cc3)c2O[C@@H]1c1ccc(O)cc1 '
               'NAME: afzelechin-(4alpha->8)-afzelechin REASON: MISSED Found 0 '
               'flavonoid core(s), need at least 2\n'
               ' * SMILES: '
               'O[C@@H]1Cc2c(O)cc(O)c([C@@H]3[C@@H](O)[C@H](Oc4cc(O)cc(O)c34)c3ccc(O)c(O)c3)c2O[C@@H]1c1ccc(O)c(O)c1 '
               'NAME: procyanidin B2 REASON: MISSED Found 0 flavonoid core(s), '
               'need at least 2\n'
               ' * SMILES: '
               'O[C@H]1Cc2c(O[C@@H]1c1ccc(O)c(O)c1)cc(O)c([C@@H]1[C@@H](O)[C@H](Oc3cc(O)cc(O)c13)c1cc(O)c(O)c(O)c1)c2O '
               'NAME: (-)-epigallocatechin-(4beta->6)-(+)-catechin REASON: '
               'MISSED Found 0 flavonoid core(s), need at least 2\n'
               ' * SMILES: '
               'O[C@@H]1Cc2c(O[C@@H]1c1ccc(O)c(O)c1)cc(O)c([C@@H]1[C@@H](O)[C@H](Oc3cc(O)cc(O)c13)c1ccc(O)c(O)c1)c2O '
               'NAME: procyanidin B5 REASON: MISSED Found 0 flavonoid core(s), '
               'need at least 2\n'
               ' * SMILES: '
               'COc1c(O)cc(cc1O)[C@H]1Oc2c(C[C@@H]1O)c(O)cc(O)c2[C@H]1[C@H](O)[C@H](Oc2cc(O)cc(O)c12)c1cc(O)c(OC)c(O)c1 '
               "NAME: (4alpha->8)-bis-4'-O-methylgallocatechin REASON: MISSED "
               'Found 0 flavonoid core(s), need at least 2\n'
               ' * SMILES: '
               'O[C@H]1Cc2c(O)cc(O)c([C@@H]3[C@@H](OC(=O)c4cc(O)c(O)c(O)c4)[C@H](Oc4cc(O)cc(O)c34)c3ccc(O)c(O)c3)c2O[C@@H]1c1ccc(O)c(O)c1 '
               'NAME: procyanidin B1 3-O-gallate REASON: MISSED Found 0 '
               'flavonoid core(s), need at least 2\n'
               ' * SMILES: '
               '[C@@]1(OC=2C=C(O)C=C(O)C2C[C@H]1O)([H])C3=CC(=O)C(=C4C(=C3)C([C@]5(OC=6C=C(O)C=C(O)C6C[C@H]5O)[H])=CC(O)=C4O)O '
               'NAME: theaflavin REASON: MISSED Found 0 flavonoid core(s), '
               'need at least 2\n'
               ' * SMILES: '
               'O[C@H]1Cc2c(O)cc(O)c([C@H]3[C@H](O)[C@H](Oc4cc(O)cc(O)c34)c3cc(O)c(O)c(O)c3)c2O[C@@H]1c1ccc(O)c(O)c1 '
               'NAME: (+)-gallocatechin-(4alpha->8)-(+)-catechin REASON: '
               'MISSED Found 0 flavonoid core(s), need at least 2\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'CNC(=O)CC[C@H](N)C(O)=O',
                                     'name': 'N(5)-methyl-L-glutamine',
                                     'reason': 'Found 0 flavonoid core(s), '
                                               'need at least 2'},
                                 {   'smiles': 'C[C@H](OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c1nc(N)[nH]c2=O)C(O)=O',
                                     'name': "L-lactyl-2-diphospho-5'-guanosine",
                                     'reason': 'Found 0 flavonoid core(s), '
                                               'need at least 2'},
                                 {   'smiles': '[H][C@@]1(CC[C@]2(C)[C@]1([H])CC[C@]1([H])[C@@]3(C)CCCC(C)(C)[C@]3([H])CC[C@@]21C)C(=C)CCC=C(C)C',
                                     'name': 'dammara-20,24-diene',
                                     'reason': 'Found 0 flavonoid core(s), '
                                               'need at least 2'},
                                 {   'smiles': 'O=C1C2=C(O)C3=C(O)C=CC=C3C(=C2C(C(=O)C)C(C1)(O)C)C=4C=5C(C(O)=C6C4C(C(=O)C)C(O)(C)CC6=O)=C(O)C=CC5',
                                     'name': 'A 39183A',
                                     'reason': 'Found 0 flavonoid core(s), '
                                               'need at least 2'},
                                 {   'smiles': 'C[C@H](O)[C@@H](O)C1=Nc2c(NC1)[nH]c(N)nc2=O',
                                     'name': 'L-threo-7,8-dihydrobiopterin',
                                     'reason': 'Found 0 flavonoid core(s), '
                                               'need at least 2'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@H]1OC)C)C)CC3=CC=C(C=C3)F',
                                     'name': 'N-[(4S,7R,8S)-5-[(4-fluorophenyl)methyl]-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'Found 0 flavonoid core(s), '
                                               'need at least 2'},
                                 {   'smiles': 'C(CCN1CCCCC1)(CCN(C(C)C)C(C)=O)(C(N)=O)C2=C(Cl)C=CC=C2',
                                     'name': 'bidisomide',
                                     'reason': 'Found 0 flavonoid core(s), '
                                               'need at least 2'},
                                 {   'smiles': 'OCCCCCCCCCCC#CC=C',
                                     'name': '13-Tetradece-11-yn-1-ol',
                                     'reason': 'Found 0 flavonoid core(s), '
                                               'need at least 2'},
                                 {   'smiles': 'C[C@@H]1CCCCO[C@@H]([C@@H](CN(C(=O)C2=C(O1)C=CC(=C2)N(C)C)[C@@H](C)CO)C)CN(C)CC3=CC=C(C=C3)OC',
                                     'name': '(3R,9S,10R)-16-(dimethylamino)-12-[(2S)-1-hydroxypropan-2-yl]-9-[[(4-methoxyphenyl)methyl-methylamino]methyl]-3,10-dimethyl-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-13-one',
                                     'reason': 'Found 0 flavonoid core(s), '
                                               'need at least 2'},
                                 {   'smiles': 'O(CCC12CC3CC(C1)CC(C2)C3)C(=O)CC4=CC=C(OCC(O)CNC(C)C)C=C4',
                                     'name': 'Adaprolol',
                                     'reason': 'Found 0 flavonoid core(s), '
                                               'need at least 2'}],
    'sample_false_negatives': [   {   'smiles': 'COc1cc2OC(CC(=O)c2c(O)c1Oc1ccc(cc1)-c1cc(=O)c2c(O)c(C)c(O)cc2o1)c1ccc(O)cc1',
                                      'name': 'taiwanhomoflavone-B',
                                      'reason': 'Flavonoid cores are not '
                                                'connected via a single atom '
                                                'or bond'},
                                  {   'smiles': 'O[C@H]1Cc2c(O)cc(O)c([C@H]3[C@H](O)[C@H](Oc4cc(O)cc(O)c34)c3ccc(O)cc3)c2O[C@@H]1c1ccc(O)cc1',
                                      'name': 'afzelechin-(4alpha->8)-afzelechin',
                                      'reason': 'Flavonoid cores are not '
                                                'connected via a single atom '
                                                'or bond'},
                                  {   'smiles': 'O[C@@H]1Cc2c(O)cc(O)c([C@@H]3[C@@H](O)[C@H](Oc4cc(O)cc(O)c34)c3ccc(O)c(O)c3)c2O[C@@H]1c1ccc(O)c(O)c1',
                                      'name': 'procyanidin B2',
                                      'reason': 'Flavonoid cores are not '
                                                'connected via a single atom '
                                                'or bond'},
                                  {   'smiles': 'O[C@H]1Cc2c(O[C@@H]1c1ccc(O)c(O)c1)cc(O)c([C@@H]1[C@@H](O)[C@H](Oc3cc(O)cc(O)c13)c1cc(O)c(O)c(O)c1)c2O',
                                      'name': '(-)-epigallocatechin-(4beta->6)-(+)-catechin',
                                      'reason': 'Flavonoid cores are not '
                                                'connected via a single atom '
                                                'or bond'},
                                  {   'smiles': 'O[C@@H]1Cc2c(O[C@@H]1c1ccc(O)c(O)c1)cc(O)c([C@@H]1[C@@H](O)[C@H](Oc3cc(O)cc(O)c13)c1ccc(O)c(O)c1)c2O',
                                      'name': 'procyanidin B5',
                                      'reason': 'Flavonoid cores are not '
                                                'connected via a single atom '
                                                'or bond'},
                                  {   'smiles': 'COc1c(O)cc(cc1O)[C@H]1Oc2c(C[C@@H]1O)c(O)cc(O)c2[C@H]1[C@H](O)[C@H](Oc2cc(O)cc(O)c12)c1cc(O)c(OC)c(O)c1',
                                      'name': "(4alpha->8)-bis-4'-O-methylgallocatechin",
                                      'reason': 'Flavonoid cores are not '
                                                'connected via a single atom '
                                                'or bond'},
                                  {   'smiles': 'O[C@H]1Cc2c(O)cc(O)c([C@@H]3[C@@H](OC(=O)c4cc(O)c(O)c(O)c4)[C@H](Oc4cc(O)cc(O)c34)c3ccc(O)c(O)c3)c2O[C@@H]1c1ccc(O)c(O)c1',
                                      'name': 'procyanidin B1 3-O-gallate',
                                      'reason': 'Flavonoid cores are not '
                                                'connected via a single atom '
                                                'or bond'},
                                  {   'smiles': '[C@@]1(OC=2C=C(O)C=C(O)C2C[C@H]1O)([H])C3=CC(=O)C(=C4C(=C3)C([C@]5(OC=6C=C(O)C=C(O)C6C[C@H]5O)[H])=CC(O)=C4O)O',
                                      'name': 'theaflavin',
                                      'reason': 'Found 1 flavonoid core(s), '
                                                'need at least 2'},
                                  {   'smiles': 'O[C@H]1Cc2c(O)cc(O)c([C@H]3[C@H](O)[C@H](Oc4cc(O)cc(O)c34)c3cc(O)c(O)c(O)c3)c2O[C@@H]1c1ccc(O)c(O)c1',
                                      'name': '(+)-gallocatechin-(4alpha->8)-(+)-catechin',
                                      'reason': 'Flavonoid cores are not '
                                                'connected via a single atom '
                                                'or bond'},
                                  {   'smiles': 'O[C@@H]1Cc2c(O)cc(O)c([C@H]3[C@@H](O)[C@H](Oc4cc(O)cc(O)c34)c3ccc(O)c(O)c3)c2O[C@@H]1c1ccc(O)c(O)c1',
                                      'name': '(-)-epicatechin-(4alpha->8)-(-)-epicatechin',
                                      'reason': 'Flavonoid cores are not '
                                                'connected via a single atom '
                                                'or bond'}],
    'attempt': 1,
    'reasoning': 'Upon analyzing the previous program and the outcomes, it '
                 'appears that the main issue lies in the way the flavonoid '
                 'core is defined using the SMARTS pattern. The SMARTS pattern '
                 "`'c1cc2oc(=O)cc2cc1'` represents a very specific "
                 'chromone-like core, which may not capture the diversity of '
                 'flavonoid structures present in biflavonoids.\n'
                 '\n'
                 'Flavonoids are a large class of polyphenolic compounds '
                 'characterized by a C6-C3-C6 skeleton, consisting of two '
                 'phenyl rings (A and B rings) and a heterocyclic ring (C '
                 'ring). Biflavonoids are dimers of flavonoids linked '
                 'together, and the connection can occur in various ways, not '
                 'limited to a single bond between chromone cores.\n'
                 '\n'
                 'To improve the program, we can:\n'
                 '\n'
                 '1. **Define a more general SMARTS pattern for flavonoid '
                 'units** that captures the various subclasses (e.g., '
                 'flavones, flavanones, flavonols, flavan-3-ols).\n'
                 '2. **Ensure that two flavonoid units are present** in the '
                 'molecule.\n'
                 '3. **Check for a connection between the two flavonoid '
                 'units** via a single atom or bond.\n'
                 '\n'
                 'Here is the updated code implementing these improvements:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 5,
    'num_false_positives': 20,
    'num_true_negatives': 142236,
    'num_false_negatives': 39,
    'num_negatives': None,
    'precision': 0.2,
    'recall': 0.11363636363636363,
    'f1': 0.14492753623188406,
    'accuracy': 0.9995853829936754,
    'negative_predictive_value': 0.9997258829731155}