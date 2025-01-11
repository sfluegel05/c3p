"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
from rdkit import Chem

def is_L_alpha_amino_acid(smiles: str) -> (bool, str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string. 
    L-alpha-amino acids have an alpha carbon with L-configuration attached to an amino group, 
    a carboxyl group, and a side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for L-alpha-amino acid pattern: C attached to N, C, and a carboxylic group
    l_alpha_amino_pattern = Chem.MolFromSmarts("[C@@H](N)(C(O)=O)")
    if not mol.HasSubstructMatch(l_alpha_amino_pattern):
        return False, "Does not match the L-alpha-amino acid pattern"
        
    return True, "Matches the L-alpha-amino acid pattern with L configuration at alpha carbon"

# Examples for testing
smiles_list = [
    "[H][C@@]1(N\C([C@@H](CCC(O)=O)[C@@H]1=)=CC1=CC(=O)c2c(C)c(Cc3[nH]c(CC4NC(=O)C(C)=C4C=C)c(C)c3CC)[nH]c12)C(O)=O",
    "C1(=CNC2=C1C=CC=C2)C([C@@H](C(=O)O)N)O",  # 3-hydroxy-L-tryptophan
    "O=C(O)[C@@H](N)CCCCCSC",  # L-trihomomethionine
    # Add other example SMILES strings here...
]

for smiles in smiles_list:
    result, reason = is_L_alpha_amino_acid(smiles)
    print(f"SMILES: {smiles}\nResult: {result}, Reason: {reason}\n")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:15705',
                          'name': 'L-alpha-amino acid',
                          'definition': 'Any alpha-amino acid having '
                                        'L-configuration at the alpha-carbon.',
                          'parents': ['CHEBI:33704'],
                          'xrefs': ['KEGG:C00151'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O=C1O[C@@H](C(=O)[C@@]2([C@]1(C(=C)C[C@H]3[C@]4([C@@H](C(C(=O)CC4)(C)C)[C@H](C[C@]23C)O)C)C)C(=O)OC)C',
                                     'name': 'Terreustoxin L',
                                     'reason': 'Does not match the '
                                               'L-alpha-amino acid pattern'},
                                 {   'smiles': '[H][C@]12CCC[C@H](C)[C@@]1(C)C[C@@H](CC2)C(C)C',
                                     'name': 'eremophilane',
                                     'reason': 'Does not match the '
                                               'L-alpha-amino acid pattern'},
                                 {   'smiles': 'N.C[C@]12CC[C@@](C)(C[C@H]1C1=CC(=O)[C@@H]3[C@@]4(C)CC[C@H](O[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@H]5O[C@@H]5O[C@@H]([C@@H](O)[C@H](O)[C@H]5O)C(O)=O)C(O)=O)C(C)(C)[C@@H]4CC[C@@]3(C)[C@]1(C)CC2)C(O)=O',
                                     'name': 'Monoammonium glycyrrhizinate',
                                     'reason': 'Does not match the '
                                               'L-alpha-amino acid pattern'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@H]2[C@H](O[C@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O)O[C@@H]1CO)[C@@H]4O[C@@H]([C@@H](O)[C@H](O[C@@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O)CO)[C@@H]4O[C@@H]6OC[C@@H](O)[C@H](O)[C@H]6O)CO',
                                     'name': 'N-[(2R,3R,4R,5S,6R)-2-[(2R,3S,4R,5R,6R)-5-Acetamido-6-hydroxy-2-(hydroxymethyl)-4-[(2R,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-4-hydroxy-5-[(2S,3S,4S,5R,6R)-5-hydroxy-6-(hydroxymethyl)-4-[(2S,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-3-[(2S,3R,4S,5R)-3,4,5-trihydroxyoxan-2-yl]oxyoxan-2-yl]oxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'Does not match the '
                                               'L-alpha-amino acid pattern'},
                                 {   'smiles': 'CN(C)CC(=O)N[C@H]1CC[C@@H](O[C@H]1CO)CCN2C=C(N=N2)C3=CC(=CC=C3)OC',
                                     'name': '2-(dimethylamino)-N-[(2R,3S,6R)-2-(hydroxymethyl)-6-[2-[4-(3-methoxyphenyl)-1-triazolyl]ethyl]-3-oxanyl]acetamide',
                                     'reason': 'Does not match the '
                                               'L-alpha-amino acid pattern'},
                                 {   'smiles': 'CC1=CC=C(C=C1)CS(=O)(=O)C2=NC=C(C(=N2)C(=O)NC3=NN=C(S3)C(C)C)Cl',
                                     'name': '5-chloro-2-[(4-methylphenyl)methylsulfonyl]-N-(5-propan-2-yl-1,3,4-thiadiazol-2-yl)-4-pyrimidinecarboxamide',
                                     'reason': 'Does not match the '
                                               'L-alpha-amino acid pattern'},
                                 {   'smiles': 'P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)CCCCCCC/C=C\\CCCCC)COC(=O)CCCCCCC/C=C\\CCCCC)([O-])=O',
                                     'name': 'PC(15:1(9Z)/15:1(9Z))',
                                     'reason': 'Does not match the '
                                               'L-alpha-amino acid pattern'},
                                 {   'smiles': 'O=C(OC[C@]1(O)[C@@H]2C=C3C=C[C@@H]4[C@@]([C@]3(C2)CC1)(CC[C@H]([C@]4(CO)C)O)C)C',
                                     'name': 'Aphidicolin A42',
                                     'reason': 'Does not match the '
                                               'L-alpha-amino acid pattern'},
                                 {   'smiles': 'O=C1OC(=CC(=C1)O)[C@H]2C3=C(C(O)=CC=C3)C(=O)C[C@]2(O)C4=CC=CC=C4',
                                     'name': 'Wailupemycin D',
                                     'reason': 'Does not match the '
                                               'L-alpha-amino acid pattern'},
                                 {   'smiles': 'COC1=C(C=C2C(=C1)C3=CC=CC=C3O2)NC(=O)CSC4=NC5=CC=CC=C5N4',
                                     'name': '2-(1H-benzimidazol-2-ylthio)-N-(2-methoxy-3-dibenzofuranyl)acetamide',
                                     'reason': 'Does not match the '
                                               'L-alpha-amino acid pattern'}],
    'sample_false_negatives': [   {   'smiles': 'S1CCC(N)(CC1)C(O)=O',
                                      'name': '4-aminotetrahydro-2H-thiopyran-4-carboxylic '
                                              'acid',
                                      'reason': 'Does not match the '
                                                'L-alpha-amino acid pattern'},
                                  {   'smiles': 'C[C@](N)(CCC(N)=O)C(O)=O',
                                      'name': '2-methyl-L-glutamine',
                                      'reason': 'Does not match the '
                                                'L-alpha-amino acid pattern'},
                                  {   'smiles': 'O=C(O)[C@@]1(N[C@@](C(=O)O)(C)C1)CC(C)C',
                                      'name': '(-)-monascumic acid',
                                      'reason': 'Does not match the '
                                                'L-alpha-amino acid pattern'},
                                  {   'smiles': 'C[C@](N)(Cc1ccc(O)cc1)C(O)=O',
                                      'name': 'alpha-methyl-L-tyrosine',
                                      'reason': 'Does not match the '
                                                'L-alpha-amino acid pattern'},
                                  {   'smiles': 'OC(=O)[C@@H]1[C@@]2([C@]1(CC[C@@]2(N)C(O)=O)[H])[H]',
                                      'name': 'Eglumetad',
                                      'reason': 'Does not match the '
                                                'L-alpha-amino acid pattern'},
                                  {   'smiles': 'C(CC(=O)[O-])[C@@]1([H])C(=N[C@](CCC(=O)[O-])([H])C(=N1)O)O',
                                      'name': 'cyclo(Glu-Glu)',
                                      'reason': 'Does not match the '
                                                'L-alpha-amino acid pattern'},
                                  {   'smiles': 'S1C(=N[C@@H](C1)[C@H]2SC[C@@](N2C)(C(=O)O)C)C3=C(O)C=CC=C3',
                                      'name': 'Thiazostatin A',
                                      'reason': 'Does not match the '
                                                'L-alpha-amino acid pattern'},
                                  {   'smiles': 'C[C@@](CC1=CC=C(OP(O)(O)=O)C(O)=C1)(NN)C(O)=O',
                                      'name': 'foscarbidopa',
                                      'reason': 'Does not match the '
                                                'L-alpha-amino acid pattern'},
                                  {   'smiles': 'OC(=O)C1(N)CCC1',
                                      'name': '1-AMINOCYCLOBUTANE CARBOXYLIC '
                                              'ACID',
                                      'reason': 'Does not match the '
                                                'L-alpha-amino acid pattern'},
                                  {   'smiles': 'O=C(O)[C@]1(N[C@](C(=O)O)(C)C1)CC(C)C',
                                      'name': '(-)-syn-2-isobutyl-4-methylazetidine-2,4-dicarboxylic '
                                              'acid',
                                      'reason': 'Does not match the '
                                                'L-alpha-amino acid pattern'}],
    'attempt': 0,
    'reasoning': 'To classify chemical entities as L-alpha-amino acids based '
                 'on their SMILES string, we need to identify specific '
                 'structural features that define this class of molecules. An '
                 'L-alpha-amino acid has the following characteristics:\n'
                 '\n'
                 '1. It contains an alpha carbon (C-alpha) with an '
                 'L-configuration. \n'
                 '2. The alpha carbon is bonded to an amino group (-NH2), a '
                 'carboxylic acid group (-COOH), and a side chain R group. \n'
                 '3. The stereochemistry at the alpha carbon is denoted by the '
                 "'@' symbol in the SMILES notation, indicating chirality, "
                 'i.e., `C@` or `C@@`.\n'
                 '\n'
                 "With this understanding, let's write a program using the "
                 'RDKit library to check if a given SMILES string represents '
                 'an L-alpha-amino acid.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 337,
    'num_false_positives': 100,
    'num_true_negatives': 909,
    'num_false_negatives': 11,
    'num_negatives': None,
    'precision': 0.7711670480549199,
    'recall': 0.9683908045977011,
    'f1': 0.8585987261146496,
    'accuracy': 0.918201915991157,
    'negative_predictive_value': 0.9880434782608696}