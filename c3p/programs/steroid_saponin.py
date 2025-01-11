"""
Classifies: CHEBI:61655 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a steroid saponin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a steroid backbone pattern (tetracyclic core)
    steroid_pattern = Chem.MolFromSmarts("C1CC2CCC3C(C2C1)CCC4C3(CCCC4)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
        
    # Look for glycosidic linkage patterns (C-O-C linkage indicative of sugars)
    glycosidic_pattern = Chem.MolFromSmarts("[C,O]-O-[C]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if len(glycosidic_matches) == 0:
        return False, "No glycosidic linkage found"
        
    # Check for typical sugar moiety (e.g., glucose, based on OH groups and specific patterns)
    sugar_pattern = Chem.MolFromSmarts("C(O)C(O)C(CO)O")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0:
        return False, "No sugar moiety found"

    return True, "Contains steroid backbone with glycosidic linkage to sugar moieties"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61655',
                          'name': 'steroid saponin',
                          'definition': 'Any saponin derived from a '
                                        'hydroxysteroid.',
                          'parents': ['CHEBI:26605', 'CHEBI:35341'],
                          'xrefs': [   'PMID:18486659',
                                       'PMID:20346608',
                                       'PMID:20846658'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'O(C1CC=2C(C3C(C4C(C5C(N6C(C5C)CCC(C6)C)C4)(CC3)C)CC2)(CC1)C)C7OC(C(O)C(O)C7OC8OC(C(O)C(O)C8O)C)CO',
                                      'name': 'beta1-Chaconine',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': 'O=C(O)[C@]12C3=C([C@@]4([C@H](C(=C)[C@@H](O[C@H]5O[C@@H]([C@@H](OC)[C@@H]([C@H]5O)O)CO)CC4)CC3)C)CC[C@@]1([C@@H]([C@H](CCC(=C)C(C)C)C)C[C@H]2O)C',
                                      'name': 'Ascosteroside',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': 'O(C1CC=2C(C3C(C4C(C5C(N6C(C5C)CCC(C6)C)C4)(CC3)C)CC2)(CC1)C)C7OC(C(OC8OC(C(O)C(OC9OCC(O)C(O)C9O)C8OC%10OC(C(O)C(O)C%10O)CO)CO)C(O)C7O)CO',
                                      'name': 'delta5-Demissine',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': 'O=C1[C@]2(O)[C@@H]([C@]3(C[C@@H](O)[C@@H]([C@@]([C@@H]3[C@H]1O)(CO[C@@H]4O[C@@H]([C@@H](O)[C@H]([C@@H]4O)O)CO)C)O)C)CC[C@@](C2)(C=C)C',
                                      'name': 'Virescenoside Z9',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': 'O1C2C(C3(C(C4C(CC3)C5(C(=CC4)CC(OC6OC(C(O)C(O)C6O)CO)CC5)C)C2)C)C(C17[NH2+]CC(CC7)C)C',
                                      'name': 'solasodine '
                                              '3-O-beta-D-glucopyranoside',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': '[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@H](CC[C@]4(C)[C@@]3([H])CC[C@]12C)O[C@@H]1O[C@H](COC(=O)CCCCCCCCCCCCCCSSCC)[C@H](O)[C@H](O)[C@H]1O)[C@H](C)CCCC(C)C',
                                      'name': 'cholesteryl '
                                              '6-O-[15-(ethyldisulfanyl)pentadecanoyl]-beta-D-galactoside',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': 'O([C@@H]1CC=2[C@@]([C@@]3([C@]([C@]4([C@](CC3)([C@H](CC4)C(=O)C)C)[H])(CC2)[H])[H])(CC1)C)[C@@H]5OC([C@@H](O[C@@H]6OC([C@H](O)C(O)[C@H]6O)C)[C@H](O)C5O)CO',
                                      'name': 'Spongipregnoloside B',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': '[H][C@@]12CC=C3C[C@H](CC[C@]3(C)[C@@]1([H])CC[C@@]1(C)[C@@]2([H])C[C@H](O[C@@H]2OC[C@H](O)[C@H](O[C@@H]3OC[C@@H](O)[C@H](O)[C@H]3O)[C@H]2OC(C)=O)[C@]1(O)[C@@H](C)C(=O)CCC(C)C)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O',
                                      'name': '3beta-[(beta-D-glucopyranosyl)oxy]-17alpha-hydroxy-16beta-[(O-beta-D-xylopyranosyl)-(1->3)-2-O-acetyl-alpha-L-arabinopyranosyloxy]cholest-5-en-22-one',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': '[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@H](CC[C@]4(C)[C@@]3([H])CC[C@]12C)O[C@@H]1O[C@H](COC(=O)CCCCC)[C@H](O)[C@H](O)[C@H]1O)[C@H](C)CCCC(C)C',
                                      'name': 'cholesteryl '
                                              '6-O-hexanoyl-beta-D-galactoside',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': 'ClCC(O)[C@@]1(CC2=CC[C@H]3[C@]([C@H](O)[C@@H](C[C@@]3([C@H]2CC1)C)O)(COC4OC(C(=O)O)C(O)C(C4O)O)C)C',
                                      'name': 'Virescenoside Z7',
                                      'reason': 'No steroid backbone found'}],
    'attempt': 0,
    'reasoning': 'To classify a chemical entity as a steroid saponin based on '
                 'its SMILES string, we need to identify characteristics '
                 'specific to that class. Steroid saponins generally consist '
                 'of a steroid backbone with one or more sugar moieties '
                 'attached. These structures often feature:\n'
                 '\n'
                 '1. A sterol-like core: composed of a fused multi-ring '
                 'structure, commonly a tetracyclic steroid nucleus.\n'
                 '2. Glycosidic linkage(s): sugar residues attached to the '
                 'steroidal backbone, often a glycosidic bond at one or '
                 'several positions.\n'
                 '\n'
                 "Given these characteristics, here's a program that aims to "
                 'classify steroid saponins by detecting the aforementioned '
                 'features in a SMILES string using RDKit.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 68,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.40476190476190477,
    'f1': 0.576271186440678,
    'accuracy': 0.40476190476190477,
    'negative_predictive_value': 0.0}