"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: CHEBI:32876 tertiary amine
A tertiary amine is a compound formally derived from ammonia by replacing three hydrogen atoms by hydrocarbyl groups.
"""
from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine has a nitrogen atom bonded to three carbon atoms (hydrocarbyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through all nitrogen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen atom
            # Check if the nitrogen is bonded to exactly 3 carbon atoms
            carbon_neighbors = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6]
            if len(carbon_neighbors) == 3:
                return True, "Contains a nitrogen atom bonded to three carbon atoms (tertiary amine)"

    # If no such nitrogen is found, it's not a tertiary amine
    return False, "No nitrogen atom bonded to three carbon atoms found"

# Example usage:
# print(is_tertiary_amine("CCN(CC)CC"))  # Should return (True, "Contains a nitrogen atom bonded to three carbon atoms (tertiary amine)")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:32876',
                          'name': 'tertiary amine',
                          'definition': 'A compound formally derived from '
                                        'ammonia by replacing three hydrogen '
                                        'atoms by hydrocarbyl groups.',
                          'parents': ['CHEBI:32952', 'CHEBI:50996'],
                          'xrefs': ['KEGG:C02196'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O1C(CCC=2C1=CC(OC3OC(C(O)C(O)C3O)C(O)=O)=C(C2OC)C=4C(=O)C=5C(OC4)=CC(O)=C(O)C5)(C)C',
                                     'name': '6-{[6-(6,7-dihydroxy-4-oxo-4H-chromen-3-yl)-5-methoxy-2,2-dimethyl-3,4-dihydro-2H-1-benzopyran-7-yl]oxy}-3,4,5-trihydroxyoxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No nitrogen atom bonded to '
                                               'three carbon atoms found'},
                                 {   'smiles': 'O=C1O[C@@H]([C@@H](OC)C=CC=C(C[C@@H](C)[C@@H]([C@@H]([C@@H]([C@@H](C=C(C=C1OC)C)C)O)C)O)C)[C@H]([C@@H](O)[C@@H]([C@@]2(O[C@H](/C=C/C)[C@@H](C)[C@@H](C2)OC3OC(C(O)C(C3)O)C)O)C)C',
                                     'name': 'Concanamycin D',
                                     'reason': 'No nitrogen atom bonded to '
                                               'three carbon atoms found'},
                                 {   'smiles': 'O1[C@@H](O[C@@H]2[C@@H](O)[C@@H](O[C@@H]([C@@H]2O)CO)O[C@H]3[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]3CO)O)[C@H](NC(=O)C)[C@@H](O[C@@H]4O[C@@H]([C@H](O)[C@H](O)[C@H]4O)CO)[C@H](O[C@@H]5O[C@H]([C@@H](O)[C@@H](O)[C@@H]5O)C)[C@H]1CO',
                                     'name': 'N-[(3R,4R,5S,6R)-5-[(2S,3R,4S,5S,6R)-4-[(2S,3R,4R,5S,6R)-3-Acetamido-6-(hydroxymethyl)-4-[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-5-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'No nitrogen atom bonded to '
                                               'three carbon atoms found'},
                                 {   'smiles': 'O=C(O)[C@]1([C@H]2[C@@](OC=3C=C4C5=C(C=CC=C5)NC4=CC3CC2)(CC[C@@H]1O)C)C',
                                     'name': 'Oxiamycin',
                                     'reason': 'No nitrogen atom bonded to '
                                               'three carbon atoms found'},
                                 {   'smiles': 'S(=O)(C(SSCCC)CC)CCC',
                                     'name': 'Propyl 1-(propylsulfinyl)propyl '
                                             'disulfide',
                                     'reason': 'No nitrogen atom bonded to '
                                               'three carbon atoms found'},
                                 {   'smiles': 'ClC=1C(=O)[C@@H]([C@@](O)(C/C=C\\CCCCC)C1)C[C@H](OC(=O)C)[C@@H](OC(=O)C)CCCC(OC)=O',
                                     'name': 'punaglandin 6',
                                     'reason': 'No nitrogen atom bonded to '
                                               'three carbon atoms found'},
                                 {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCCCC[C@H](O)C(=O)N[C@@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[C@H](O)[C@@H](O)CCCCCCCCCCCCCC',
                                     'name': 'Ins-1-P-Cer(t18:0/2-OH-26:0)',
                                     'reason': 'No nitrogen atom bonded to '
                                               'three carbon atoms found'},
                                 {   'smiles': 'CCC(C)(C)C(=O)OC1CC(C=C2C1[C@H]([C@H](C=C2)C)CC[C@@H]3CC(CC(=O)O3)O)C',
                                     'name': '2,2-dimethylbutanoic acid '
                                             '[(7S,8S)-8-[2-[(2R)-4-hydroxy-6-oxo-2-oxanyl]ethyl]-3,7-dimethyl-1,2,3,7,8,8a-hexahydronaphthalen-1-yl] '
                                             'ester',
                                     'reason': 'No nitrogen atom bonded to '
                                               'three carbon atoms found'},
                                 {   'smiles': '[C@@H]1(C2=CC=CC=C2)[C@@H](N)C1',
                                     'name': '(1S,2R)-tranylcypromine',
                                     'reason': 'No nitrogen atom bonded to '
                                               'three carbon atoms found'},
                                 {   'smiles': 'OC(=O)/C=C/C#CC#CC#N',
                                     'name': '(2E)-7-Cyanohept-2-en-4,6-diynoic '
                                             'acid',
                                     'reason': 'No nitrogen atom bonded to '
                                               'three carbon atoms found'}],
    'sample_false_negatives': [   {   'smiles': 'O=C(/N=[N+](\\[O-])/C1=CC=C(O)C=C1)N',
                                      'name': '4-hydroxyphenyl-1-ONN-azoxyformamide',
                                      'reason': 'No nitrogen atom bonded to '
                                                'three carbon atoms found'},
                                  {   'smiles': 'O=C(O)/C=C/C=C/C1=CC=C(N=[N+]([O-])C2=CC=C(/C=C/C=C/C(=O)N[C@H](C(=O)O)CCC(=O)N)C=C2)C=C1',
                                      'name': 'Azoxymycin B',
                                      'reason': 'No nitrogen atom bonded to '
                                                'three carbon atoms found'},
                                  {   'smiles': 'O=C(NC(C(=O)O)CCC(=O)N)/C=C/C=C/C1=CC=C([N+]([O-])=NC2=CC=C(/C=C/C=C/C(=O)NC(C(=O)O)CCC(=O)N)C=C2)C=C1',
                                      'name': 'Azoxymycin A',
                                      'reason': 'No nitrogen atom bonded to '
                                                'three carbon atoms found'},
                                  {   'smiles': 'C/C(=C\\CC\\C(\\C)=C\\C[C@]12C(C=C(C([C@H]1O2)OC(C)=O)C)=O)/CCCC(C)(C)O',
                                      'name': 'yanuthone M',
                                      'reason': 'No nitrogen atom bonded to '
                                                'three carbon atoms found'}],
    'attempt': 0,
    'reasoning': 'To classify a chemical entity as a tertiary amine, we need '
                 'to identify a nitrogen atom that is bonded to three carbon '
                 'atoms (hydrocarbyl groups) and no hydrogen atoms. This can '
                 'be done by checking the connectivity of nitrogen atoms in '
                 'the molecule.\n'
                 '\n'
                 'Here’s the Python program that accomplishes this:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 72,
    'num_false_positives': 100,
    'num_true_negatives': 287,
    'num_false_negatives': 4,
    'num_negatives': None,
    'precision': 0.4186046511627907,
    'recall': 0.9473684210526315,
    'f1': 0.5806451612903226,
    'accuracy': 0.775377969762419,
    'negative_predictive_value': 0.9862542955326461}