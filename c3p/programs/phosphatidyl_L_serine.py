"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
"""
Classifies: phosphatidyl-L-serine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a phosphatidyl-L-serine based on its SMILES string.
    A phosphatidyl-L-serine is a glycerophospholipid with a phosphatidyl group esterified
    to the hydroxy group of L-serine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone (3 connected carbon atoms)
    glycerol_pattern = Chem.MolFromSmarts("[CH2][CH][CH2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for at least 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Check for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for serine moiety connected to phosphate
    # Serine pattern: N[C@H](C(=O)O)CO
    # Phosphate connected to serine hydroxy group
    serine_phosphate_pattern = Chem.MolFromSmarts("N[C@H](C(=O)O)CO[P](=O)(O)O")
    if not mol.HasSubstructMatch(serine_phosphate_pattern):
        return False, "No serine moiety connected to phosphate group found"

    # All criteria met
    return True, "Contains glycerol backbone with two fatty acid chains, phosphate group, and serine moiety"

__metadata__ = {
    'chemical_class': {
        'name': 'phosphatidyl-L-serine',
        'definition': 'A class of aminophospholipids in which a phosphatidyl group is esterified to the hydroxy group of serine.'
    }
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18303',
                          'name': 'phosphatidyl-L-serine',
                          'definition': 'A class of aminophospholipids in '
                                        'which a phosphatidyl group is '
                                        'esterified to the hydroxy group of '
                                        'serine.',
                          'parents': [   'CHEBI:52565',
                                         'CHEBI:60971',
                                         'CHEBI:84135'],
                          'xrefs': [   'DrugBank:DB00144',
                                       'HMDB:HMDB0014291',
                                       'KEGG:C02737',
                                       'MetaCyc:L-1-PHOSPHATIDYL-SERINE',
                                       'PMID:10540156',
                                       'PMID:15533308',
                                       'PMID:19687511',
                                       'PMID:23543734',
                                       'PMID:3106116',
                                       'PMID:3196084',
                                       'PMID:4153523',
                                       'PMID:8204602',
                                       'PMID:8626656',
                                       'PMID:9677350',
                                       'Patent:EP2322184',
                                       'Patent:HK1046237',
                                       'Patent:US2011098249',
                                       'Wikipedia:Phosphatidylserine'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C3=CC=C(C=C3)C(=O)N(C)C)O[C@@H]1CN(C)C(=O)NC4=CC=CC=C4F)[C@H](C)CO',
                                     'name': '4-[(4R,5S)-5-[[[(2-fluoroanilino)-oxomethyl]-methylamino]methyl]-2-[(2R)-1-hydroxypropan-2-yl]-4-methyl-1,1-dioxo-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-8-yl]-N,N-dimethylbenzamide',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'C[C@H]1[C@H](CN(C(=O)C2=C(C=CC(=C2)C#N)OC[C@H]3[C@H](CC[C@@H](O3)CCN(C1=O)C)OC)C)OC',
                                     'name': 'LSM-10564',
                                     'reason': 'Found 0 ester groups, need at '
                                               'least 2'},
                                 {   'smiles': 'CCC(C(=O)NC1=NN=C(S1)SCC)(C(F)(F)F)C(F)(F)F',
                                     'name': 'N-[5-(ethylthio)-1,3,4-thiadiazol-2-yl]-2,2-bis(trifluoromethyl)butanamide',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O=C(N[C@@H](CCC(=O)N)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC1=CC=C(O)C=C1)CC(=O)N',
                                     'name': 'Tyr-Asn-Gln',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'OC[C@H]1OC(O)[C@@H](O)[C@@H]1O',
                                     'name': 'D-arabinofuranose',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O[C@@H]([C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\\3/C[C@@H](O)CCC3=C)[H])C)[H])C)C',
                                     'name': '(5Z,7E)-(3S,22R)- '
                                             '24-nor-9,10-seco-5,7,10(19)-cholatriene-3,22-diol',
                                     'reason': 'Found 0 ester groups, need at '
                                               'least 2'},
                                 {   'smiles': 'O1C(C(CC1=O)C)CCCCCC',
                                     'name': '5-Hexyldihydro-4-methyl-2(3H)-furanone',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O=C1O[C@]23C(=O)N[C@H]([C@@H]2[C@H](C)C(=C[C@@H]3C=CC[C@@H]([C@H]([C@H](C=C1)C)O)C)C)CC4=CC=CC=C4',
                                     'name': '7-deoxy-cytochalasin Z7',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': '[C@]123[C@H](O1)[C@](C(C=C2[C@]4([C@](CC3)(C([C@@H](C(C4)=O)O)(C)C)[H])C)=O)(C=C)C',
                                     'name': 'yucalexin P8',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': '[H][C@]12CC[C@@]3([H])[C@]4([H])CC[C@]([H])([C@H](C)CCC(O)=O)[C@@]4(C)[C@@H](O)C[C@]3([H])[C@@]1(C)CC[C@@H](O)C2',
                                     'name': 'deoxycholic acid',
                                     'reason': 'Found 0 ester groups, need at '
                                               'least 2'}],
    'sample_false_negatives': [   {   'smiles': 'P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCC/C=C\\C/C=C\\CCCC)(O)=O',
                                      'name': 'PI(17:2(9Z,12Z)/20:5(5Z,8Z,11Z,14Z,17Z))',
                                      'reason': 'No serine moiety connected to '
                                                'phosphate group found'}],
    'attempt': 0,
    'reasoning': 'To determine whether a molecule is a phosphatidyl-L-serine, '
                 'we need to check for the following features:\n'
                 '\n'
                 '1. **Glycerol Backbone**: A phosphatidyl-L-serine has a '
                 'glycerol backbone consisting of three carbon atoms connected '
                 'sequentially.\n'
                 '2. **Ester Bonds**: There should be two fatty acid chains '
                 'attached to the glycerol backbone via ester bonds at the '
                 'sn-1 and sn-2 positions.\n'
                 '3. **Phosphate Group**: A phosphate group attached to the '
                 'sn-3 position of the glycerol backbone.\n'
                 '4. **Serine Moiety**: The phosphate group is esterified to '
                 'the hydroxy group of L-serine, which includes an amino group '
                 'and a carboxylic acid.\n'
                 '\n'
                 'The program will:\n'
                 '\n'
                 '- Parse the SMILES string.\n'
                 '- Check for the glycerol backbone.\n'
                 '- Verify the presence of at least two ester groups.\n'
                 '- Confirm the presence of a phosphate group.\n'
                 '- Check for the serine moiety connected to the phosphate '
                 'group.\n'
                 '\n'
                 "Here's the Python code implementing the above logic:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 317,
    'num_false_positives': 29,
    'num_true_negatives': 141953,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.9161849710982659,
    'recall': 0.9968553459119497,
    'f1': 0.9548192771084337,
    'accuracy': 0.9997891777933943,
    'negative_predictive_value': 0.9999929554644462}