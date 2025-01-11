"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: secondary alcohol
"""
from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.
    
    A secondary alcohol is a compound in which a hydroxy group, -OH, is attached 
    to a saturated carbon atom which has two other carbon atoms attached to it.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alcohol, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a secondary alcohol
    # [C;X4;H1;D3] - sp3 carbon with one hydrogen and degree 3
    # Attached to [O;H] - hydroxyl group
    secondary_alcohol_pattern = Chem.MolFromSmarts('[C;X4;H1;D3]([O;H])')

    # Search for matches of the secondary alcohol pattern
    matches = mol.GetSubstructMatches(secondary_alcohol_pattern)
    if matches:
        return True, "Contains a secondary alcohol group"
    else:
        return False, "No secondary alcohol group found"

__metadata__ = {
    'chemical_class': {
        'name': 'secondary alcohol',
        'definition': 'A secondary alcohol is a compound in which a hydroxy group, -OH, '
                      'is attached to a saturated carbon atom which has two other carbon '
                      'atoms attached to it.'
    }
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35681',
                          'name': 'secondary alcohol',
                          'definition': 'A secondary alcohol is a compound in '
                                        'which a hydroxy group, -OH, is '
                                        'attached to a saturated carbon atom '
                                        'which has two other carbon atoms '
                                        'attached to it.',
                          'parents': ['CHEBI:30879'],
                          'xrefs': ['KEGG:C00432', 'KEGG:C01612'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'CNC(=O)CC[C@H](N)C(O)=O',
                                     'name': 'N(5)-methyl-L-glutamine',
                                     'reason': 'No secondary alcohol group '
                                               'found'},
                                 {   'smiles': '[H][C@@]1(CC[C@]2(C)[C@]1([H])CC[C@]1([H])[C@@]3(C)CCCC(C)(C)[C@]3([H])CC[C@@]21C)C(=C)CCC=C(C)C',
                                     'name': 'dammara-20,24-diene',
                                     'reason': 'No secondary alcohol group '
                                               'found'},
                                 {   'smiles': 'O=C1C2=C(O)C3=C(O)C=CC=C3C(=C2C(C(=O)C)C(C1)(O)C)C=4C=5C(C(O)=C6C4C(C(=O)C)C(O)(C)CC6=O)=C(O)C=CC5',
                                     'name': 'A 39183A',
                                     'reason': 'No secondary alcohol group '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@H]1OC)C)C)CC3=CC=C(C=C3)F',
                                     'name': 'N-[(4S,7R,8S)-5-[(4-fluorophenyl)methyl]-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'No secondary alcohol group '
                                               'found'},
                                 {   'smiles': 'C(CCN1CCCCC1)(CCN(C(C)C)C(C)=O)(C(N)=O)C2=C(Cl)C=CC=C2',
                                     'name': 'bidisomide',
                                     'reason': 'No secondary alcohol group '
                                               'found'},
                                 {   'smiles': 'OCCCCCCCCCCC#CC=C',
                                     'name': '13-Tetradece-11-yn-1-ol',
                                     'reason': 'No secondary alcohol group '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1CCCCO[C@@H]([C@@H](CN(C(=O)C2=C(O1)C=CC(=C2)N(C)C)[C@@H](C)CO)C)CN(C)CC3=CC=C(C=C3)OC',
                                     'name': '(3R,9S,10R)-16-(dimethylamino)-12-[(2S)-1-hydroxypropan-2-yl]-9-[[(4-methoxyphenyl)methyl-methylamino]methyl]-3,10-dimethyl-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-13-one',
                                     'reason': 'No secondary alcohol group '
                                               'found'},
                                 {   'smiles': 'CC1=CC=CC=C1C2=CC=C(C=C2)C3[C@H]4CNC[C@@H]3N4C(=O)NC5=CC=C(C=C5)Cl',
                                     'name': '(1S,5R)-N-(4-chlorophenyl)-7-[4-(2-methylphenyl)phenyl]-3,6-diazabicyclo[3.1.1]heptane-6-carboxamide',
                                     'reason': 'No secondary alcohol group '
                                               'found'},
                                 {   'smiles': 'OC(=O)c1cccc(Oc2ccccc2)c1',
                                     'name': '3-phenoxybenzoic acid',
                                     'reason': 'No secondary alcohol group '
                                               'found'},
                                 {   'smiles': 'COC(=O)C1=C(C)NC(C)=C(C1c1ccccc1[N+]([O-])=O)C(=O)OCC(C)C',
                                     'name': 'methyl 2-methylpropyl '
                                             '2,6-dimethyl-4-(2-nitrophenyl)-1,4-dihydropyridine-3,5-dicarboxylate',
                                     'reason': 'No secondary alcohol group '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'O=C1C=C(C#N)C[C@]1(O)C=C',
                                      'name': 'Homothallin II',
                                      'reason': 'No secondary alcohol group '
                                                'found'},
                                  {   'smiles': 'CC(O)(Cc1ccc(O)cc1)c1ccc(O)cc1',
                                      'name': '1,2-bis(4-hydroxyphenyl)propan-2-ol',
                                      'reason': 'No secondary alcohol group '
                                                'found'},
                                  {   'smiles': 'OC(C)(C)C(=O)C',
                                      'name': '3-Hydroxy-3-methylbutan-2-one',
                                      'reason': 'No secondary alcohol group '
                                                'found'},
                                  {   'smiles': 'OC1(CCCC(C1=O)(C)C)C',
                                      'name': '2-Hydroxy-2,6,6-trimethylcyclohexanone',
                                      'reason': 'No secondary alcohol group '
                                                'found'},
                                  {   'smiles': 'BrC=1C(=O)[C@](O)(C=C)CC1N',
                                      'name': 'Bromomyrothenone B',
                                      'reason': 'No secondary alcohol group '
                                                'found'},
                                  {   'smiles': '[C@@]123[C@@]4([C@]5([C@@](CC1=O)(C(CN(CC2)C)=CCO[C@]5(CC(N4C6=C3C=CC=C6)=O)[H])O)[H])[H]',
                                      'name': '14-hydroxyicajine',
                                      'reason': 'No secondary alcohol group '
                                                'found'},
                                  {   'smiles': 'O=C1/C(=C(\\O)/C=C/C=C/C)/[C@H]2[C@@](O)(C)C([C@@]1([C@H]3C[C@@](O)(CO)O[C@@H]23)C)=O',
                                      'name': 'Sorbicillfuran A',
                                      'reason': 'No secondary alcohol group '
                                                'found'},
                                  {   'smiles': 'O=C1C(=CC([C@@]1(O)C[C@H]2C(=C)CC[C@@H]3[C@@]2(CCCC3(C)C)C)=O)C',
                                      'name': 'Penicilliumin B',
                                      'reason': 'No secondary alcohol group '
                                                'found'},
                                  {   'smiles': 'O=C1C=C[C@H]([C@]1(O)C(=O)CCCCCCCCCCCC)OC(=O)C',
                                      'name': '4-O-acetyl hygrophorone C12',
                                      'reason': 'No secondary alcohol group '
                                                'found'},
                                  {   'smiles': 'O=C1C(O)=C(OC)[C@@]([C@@]1(O)C)(O)C(=O)OC',
                                      'name': '2,4,5-trihydroxy-3-methoxy-4-methoxycarbonyl-5-methyl-2-cyclopenten-1-one',
                                      'reason': 'No secondary alcohol group '
                                                'found'}],
    'attempt': 0,
    'reasoning': 'To determine if a molecule is a **secondary alcohol**, we '
                 'need to check for the presence of at least one hydroxyl '
                 'group (-OH) attached to a saturated carbon atom that has '
                 'exactly two other carbon atoms attached to it.\n'
                 '\n'
                 "Here are the steps we'll follow:\n"
                 '\n'
                 '1. **Parse the SMILES string** to create an RDKit molecule '
                 'object.\n'
                 '2. **Validate the molecule** to ensure it was parsed '
                 'correctly.\n'
                 '3. **Define a SMARTS pattern** that matches secondary '
                 'alcohols:\n'
                 '   - A carbon atom (`C`) that is sp³ hybridized (`X4`), has '
                 'one hydrogen (`H1`), and is connected to exactly two other '
                 'carbons.\n'
                 '   - This carbon is connected to a hydroxyl group '
                 '(`[O;H]`).\n'
                 '4. **Search for matches** of the secondary alcohol pattern '
                 'in the molecule.\n'
                 '5. **Return the classification result** along with the '
                 'reason.\n'
                 '\n'
                 'Now, here is the code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 1811,
    'num_false_positives': 100,
    'num_true_negatives': 225,
    'num_false_negatives': 32,
    'num_negatives': None,
    'precision': 0.9476713762428048,
    'recall': 0.9826370048833424,
    'f1': 0.9648375066595631,
    'accuracy': 0.9391143911439115,
    'negative_predictive_value': 0.8754863813229572}