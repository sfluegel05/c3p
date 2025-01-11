"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: anthocyanidin cation
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    Anthocyanidin cations are oxygenated derivatives of flavylium (2-phenylchromenylium).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for positive charge
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != 1:
        return False, "Must have a total charge of +1"
        
    # Look for flavylium core structure (2-phenylchromenylium)
    # [O+]=C1C=CC=CC=C1 connected to phenyl ring
    flavylium_pattern = Chem.MolFromSmarts('[O+]=C1c2ccccc2C=Cc3ccccc13')
    if not mol.HasSubstructMatch(flavylium_pattern):
        return False, "No flavylium core structure found"
        
    # Must have oxygen atoms (typically hydroxyls) on the rings
    oxygen_pattern = Chem.MolFromSmarts('cO')
    oxygen_matches = mol.GetSubstructMatches(oxygen_pattern)
    if len(oxygen_matches) < 2:
        return False, "Insufficient hydroxyl/oxygen substituents"
        
    # Check for characteristic benzopyran structure
    benzopyran_pattern = Chem.MolFromSmarts('c1c2c(cc1)OC=CC2')
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No benzopyran structure found"
        
    # Count number of rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Must have at least 3 rings"
        
    # Look for characteristic hydroxyl pattern on A and B rings
    hydroxyl_pattern = Chem.MolFromSmarts('c(O)c(O)')
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing characteristic hydroxyl pattern"
        
    # Additional checks for substituents commonly found in anthocyanidins
    # (optional methoxy groups)
    methoxy_pattern = Chem.MolFromSmarts('cOC')
    has_methoxy = mol.HasSubstructMatch(methoxy_pattern)
    
    # Success message varies based on presence of methoxy groups
    base_message = "Contains flavylium core with characteristic hydroxyl pattern"
    if has_methoxy:
        message = base_message + " and methoxy substituents"
    else:
        message = base_message
        
    return True, message


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16366',
                          'name': 'anthocyanidin cation',
                          'definition': 'Any organic cation that is an aglycon '
                                        'of anthocyanin cation; they are '
                                        'oxygenated derivatives of flavylium '
                                        '(2-phenylchromenylium).',
                          'parents': ['CHEBI:25697', 'CHEBI:47916'],
                          'xrefs': [   'HMDB:HMDB0031460',
                                       'KEGG:C02003',
                                       'MetaCyc:Anthocyanidins',
                                       'Wikipedia:Anthocyanidins'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'O1C([C@@H](O)C(O)C(O)[C@@H]1OC=2C(=[O+]C=3C(C2)=C(O)C=C(O)C3)C4=CC(O)=C(O)C=C4)COC(=O)C(O)=O',
                                      'name': 'Ophrysanin',
                                      'reason': 'No flavylium core structure '
                                                'found'},
                                  {   'smiles': 'OC[C@H]1O[C@@H](Oc2cc3c(O)cc(O)cc3[o+]c2-c2ccc(O)c(O)c2)[C@H](O[C@@H]2OC[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H]1O',
                                      'name': 'cyanidin '
                                              '3-O-[beta-D-xylosyl-(1->2)-beta-D-galactoside]',
                                      'reason': 'No flavylium core structure '
                                                'found'},
                                  {   'smiles': 'OC[C@H]1O[C@@H](O[C@H]2[C@H](Oc3cc4c(O)cc(O[C@@H]5O[C@H](COC(=O)\\C=C\\c6ccc(O)c(O)c6)[C@@H](O)[C@H](O)[C@H]5O)cc4[o+]c3-c3ccc(O)c(O[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)C(O)=O)c3)O[C@H](COC(=O)CC(=O)OC(C(O)C(O)=O)C(O)=O)[C@H](O)[C@@H]2O)[C@H](OC(=O)\\C=C\\c2ccc(O)c(O)c2)[C@@H](O)[C@@H]1O',
                                      'name': 'Anemone purple anthocyanin 1',
                                      'reason': 'No flavylium core structure '
                                                'found'},
                                  {   'smiles': 'O(C1C(O)C(O)C(OC1OC=2C=C3C(OC4OC(C(O)C(O)C4O)CO)=CC(O)=CC3=[O+]C2C5=CC=C(O)C=C5)CO)C6OC(C(O)C(O)C6O)CO',
                                      'name': 'Pelargonidin 3-sophoroside '
                                              '5-glucoside',
                                      'reason': 'No flavylium core structure '
                                                'found'},
                                  {   'smiles': 'OC[C@H]1O[C@@H](Oc2cc3c(O)cc(O)cc3[o+]c2-c2ccc(O)cc2)[C@H](O[C@@H]2OC[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@@H]1O',
                                      'name': 'pelargonidin '
                                              '3-O-beta-D-sambubioside',
                                      'reason': 'No flavylium core structure '
                                                'found'},
                                  {   'smiles': 'Oc1cc(O)c2cc(O)c([o+]c2c1)-c1ccc(O)c(O)c1',
                                      'name': 'cyanidin cation',
                                      'reason': 'No flavylium core structure '
                                                'found'},
                                  {   'smiles': 'O1C([C@@H](O)[C@H](O)C(O)[C@@H]1OCC2O[C@@H](OC=3C(=[O+]C=4C(C3)=C(O)C=C(O)C4)C=5C=C(OC)C(O)=C(O)C5)C(O)C(O)[C@@H]2O)CO[C@@H]6OC([C@@H](O)[C@H](O)C6O)CO',
                                      'name': 'Petunidin 3-gentiotrioside',
                                      'reason': 'No flavylium core structure '
                                                'found'},
                                  {   'smiles': 'O([C@@H]1[C@@H](O)[C@H](O)[C@H](O[C@@H]1OC=2C=C3C(O[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)CO)=CC(O)=CC3=[O+]C2C5=CC(OC)=C(O)C=C5)COC(=O)C6=CC=C(O)C=C6)[C@@H]7O[C@@H]([C@@H](O)[C@H](O)[C@H]7O)CO',
                                      'name': 'Peonidin 3-p-hydroxybenzoyl '
                                              'sophoroside-5-glucoside',
                                      'reason': 'No flavylium core structure '
                                                'found'},
                                  {   'smiles': 'O1C([C@@H](O)C(O)C(O)[C@@H]1OC=2C(=[O+]C=3C(C2)=C(O)C=C(O)C3)C4=CC(OC)=C(O)C=C4)COC(=O)/C=C/C5=CC=C(O)C=C5',
                                      'name': 'Peonidin '
                                              "3-(6''-p-coumarylglucoside)",
                                      'reason': 'No flavylium core structure '
                                                'found'},
                                  {   'smiles': 'COc1cc(ccc1O)-c1[o+]c2cc(O)cc(O)c2cc1O',
                                      'name': 'peonidin',
                                      'reason': 'No flavylium core structure '
                                                'found'}],
    'attempt': 0,
    'reasoning': "I'll create a program to classify anthocyanidin cations. "
                 'From analyzing the examples and the definition, I can '
                 'identify key structural features:\n'
                 '\n'
                 '1. They contain a flavylium (2-phenylchromenylium) core '
                 'structure\n'
                 '2. They have a positive charge\n'
                 '3. They contain multiple hydroxyl groups\n'
                 '4. They may have additional substituents like methoxy '
                 'groups\n'
                 '5. They must be the aglycone part of anthocyanins (though '
                 'many examples show glycosylated forms)\n'
                 '\n'
                 "Here's the program:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.0,
    'negative_predictive_value': 0.0}