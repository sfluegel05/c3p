"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    A cyclic fatty acid contains a fatty acid chain with a carboxylic acid group and a ring structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for a ring structure
    if not mol.GetRingInfo().NumRings():
        return False, "No ring structure found"
    
    return True, "Contains a carboxylic acid group and a ring structure, classifying as a cyclic fatty acid"

# Example usage:
# result, reason = is_cyclic_fatty_acid("C1(C(C/C=C\CCCCCO)O1)CCCCCCCC(=O)O")
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59238',
                          'name': 'cyclic fatty acid',
                          'definition': 'Any fatty acid containing anywhere in '
                                        'its structure a ring of atoms.',
                          'parents': ['CHEBI:35366'],
                          'xrefs': ['PMID:17113094', 'PMID:9300789'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'FC(F)(F)c1nc2c(Cl)c(Cl)c(Cl)c(Cl)c2[nH]1',
                                     'name': '4,5,6,7-Tetrachloro-2-trifluoromethylbenzimidazole',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'Cl.OC[C@H]1O[C@@H]2[C@@H](Oc3nc(=N)ccn23)[C@@H]1O',
                                     'name': 'ancitabine hydrochloride',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C3=CC=CC(=C3)C#N)O[C@H]1CN(C)C(=O)C4=CC=NC=C4)[C@H](C)CO',
                                     'name': 'N-[[(4R,5R)-8-(3-cyanophenyl)-2-[(2R)-1-hydroxypropan-2-yl]-4-methyl-1,1-dioxo-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-5-yl]methyl]-N-methyl-4-pyridinecarboxamide',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'OC1=C(NC(=O)NC2=CC=CC=C2Br)C=CC(=C1)[N+]([O-])=O',
                                     'name': 'SB225002',
                                     'reason': 'No carboxylic acid group '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'O=C(C=1OC([C@H](O)[C@@H](O)CCCCC)=CC1)C',
                                      'name': 'Armillariol C',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'O1C(CCCC)=CC=C1CCC2=CC(OC)=C(O)C=C2',
                                      'name': '2-Butyl-5-[2-(4-hydroxy-3-methoxyphenyl)ethyl]furan',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'O1C(=CC=C1C(OC)CCC=C)C(O)C(O)C',
                                      'name': '1-(5-(1-methoxypent-4-en-1-yl)furan-2-yl)propane-1,2-diol',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'O=C(C=1OC(C(O)C(O)C)=CC1)CCCC',
                                      'name': '1-(5-(1,2-dihydroxypropyl)-furan-2-yl)-pentan-1-one',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'O1C(=CC=C1CC(O)CC)CC(O)C(OC)C',
                                      'name': '1-(5-(2-hydroxybutyl)furan-2-yl)-3-methoxybutan-2-OL',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'S(/C=C\\C#CCC=1OC(=CC1)/C=C\\C=O)C',
                                      'name': '(2E,11Z)-5-[5-(Methylthio)-4-penten-2-ynyl]-2-furanacrolein',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'O=C(C=1OC(C(O)CCCC)=CC1)C(O)C',
                                      'name': '2-hydroxy-1-(5-(1-hydroxypentyl)-furan-2-yl)-propan-1-one',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'O=C1C2=C(OC(=C2C[C@@H]([C@H]1O)O)C(=O)C(=O)C3=CC=CC=C3)/C=C(/OC)\\C',
                                      'name': 'Neoenterocin B',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'O=C(C=1OC(CCCC(=O)C)=CC1)[C@@H](O)C',
                                      'name': 'Irpexlacte D',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'O1C(=CC=C1C(OC)CCCC)C(O)C(O)C',
                                      'name': '1-(5-(1-methoxypentyl)furan-2-yl)propane-1,2-diol',
                                      'reason': 'No carboxylic acid group '
                                                'found'}],
    'attempt': 0,
    'reasoning': 'To classify chemical entities as cyclic fatty acids based on '
                 'their SMILES strings, we need a program that identifies the '
                 'presence of a fatty acid chain and at least one ring '
                 'structure within the molecule. A fatty acid typically '
                 'contains a long hydrocarbon chain and a carboxylic acid '
                 'group. For cyclic fatty acids, the molecule also needs to '
                 'contain a ring structure.\n'
                 '\n'
                 '1. **Import RDKit:** We will use the RDKit library, which '
                 'includes methods for molecular structure parsing and '
                 'substructure searching.\n'
                 '\n'
                 '2. **Parse the SMILES string:** Convert the SMILES '
                 'representation into a molecule object.\n'
                 '\n'
                 '3. **Check for the carboxylic acid group:** Identify the '
                 'presence of the carboxylic acid group (`-C(=O)O`) at the end '
                 'of the alkyl chain.\n'
                 '\n'
                 '4. **Identify the presence of a ring:** Use RDKit to find '
                 'ring structures within the molecule.\n'
                 '\n'
                 '5. **Return the results:** If both a ring and a fatty acid '
                 'structure (with a carboxylic acid group) are found, classify '
                 'the molecule as a cyclic fatty acid.\n'
                 '\n'
                 "Here's the implementation of the classification function:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 202,
    'num_false_positives': 100,
    'num_true_negatives': 274,
    'num_false_negatives': 12,
    'num_negatives': None,
    'precision': 0.6688741721854304,
    'recall': 0.9439252336448598,
    'f1': 0.7829457364341085,
    'accuracy': 0.8095238095238095,
    'negative_predictive_value': 0.958041958041958}