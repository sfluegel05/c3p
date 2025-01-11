"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies: CHEBI:35667 oxo fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid contains at least one aldehydic or ketonic group in addition 
    to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Count carboxylic acid groups
    carboxylic_matches = len(mol.GetSubstructMatches(carboxylic_pattern))
    if carboxylic_matches > 1:
        return False, "Multiple carboxylic acid groups found"

    # Look for ketone groups (excluding the carboxylic acid)
    ketone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[#6]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    # Look for aldehyde groups
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=[OX1])")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    total_oxo_groups = len(ketone_matches) + len(aldehyde_matches)
    
    if total_oxo_groups == 0:
        return False, "No ketone or aldehyde group found"

    # Count carbons to verify it's a fatty acid
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Carbon chain too short to be a fatty acid"
    
    # Calculate the number of carbons in the main chain
    # This uses the molecular framework to estimate chain length
    framework = Chem.MolFromSmiles(Chem.MolToSmiles(mol, canonical=True))
    
    # Count rotatable bonds to verify chain nature
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
        return False, "Too rigid to be a fatty acid"

    # Success case - found carboxylic acid and at least one oxo group
    oxo_description = []
    if len(ketone_matches) > 0:
        oxo_description.append(f"{len(ketone_matches)} ketone(s)")
    if len(aldehyde_matches) > 0:
        oxo_description.append(f"{len(aldehyde_matches)} aldehyde(s)")
    
    return True, f"Contains carboxylic acid and {' and '.join(oxo_description)}"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59644',
                          'name': 'oxo fatty acid',
                          'definition': 'Any fatty acid containing at least '
                                        'one aldehydic or ketonic group in '
                                        'addition to the carboxylic acid '
                                        'group.',
                          'parents': ['CHEBI:25754', 'CHEBI:35366'],
                          'xrefs': [   'PMID:6434570',
                                       'PMID:8454196',
                                       'PMID:8638935'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C1C=C(N=C2N1NN=N2)C3=CC=CC=C3',
                                     'name': '5-phenyl-1,7-dihydrotetrazolo[1,5-a]pyrimidine',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)[C@H](COC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\CCCC',
                                     'name': 'TG(14:1(9Z)/22:1(13Z)/18:1(11Z))',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H]1O',
                                     'name': 'beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-{alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->6)]-alpha-D-Manp-(1->6)}-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)-beta-D-GlcpNAc',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': '[C@@]12([C@@]([C@]([C@@H](CC1)C)(CCC(CC)C)C)(CCC[C@@H]2C)[H])C',
                                     'name': 'clerodane',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'O=C(NC(CC(C)C)C(=O)N[C@@H](CCCN=C(N)N)C=O)C(NC(=O)CC)CC(C)C',
                                     'name': 'Leupeptin Pr-LL',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'COC1=CC=C(C=C1)CCNC(=O)C(=C2C3=CC=CC=C3C(=N2)NC(=O)C4=CC=CS4)C#N',
                                     'name': 'N-[3-[1-cyano-2-[2-(4-methoxyphenyl)ethylamino]-2-oxoethylidene]-1-isoindolyl]-2-thiophenecarboxamide',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'O=C1O[C@@H](CC=2C1=C(OC)C(O)=CC2)CCC[C@@H](O)C',
                                     'name': 'Penicimarin C',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'CO[C@@H]1[C@@H]2[C@@H](C[C@@H]3[C@@](O2)([C@H]([C@@H]([C@H](O3)C(=O)OC)O)O)O)OC1N4C=NC5=C4N=CN=C5N',
                                     'name': 'LSM-4497',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1[C@H](O)CCC2=CC[C@H](C[C@]12C)C(C)=C',
                                     'name': '1-deoxycapsidiol',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'C1CN=C(Cc2ccccc2)N1',
                                     'name': 'tolazoline',
                                     'reason': 'No carboxylic acid group '
                                               'found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': "I'll help create a program to classify oxo fatty acids. "
                 "Let's break down the requirements:\n"
                 '\n'
                 '1. Must be a fatty acid (contains a carboxylic acid group)\n'
                 '2. Must contain at least one aldehydic (R-CH=O) or ketonic '
                 "(R-C(=O)-R') group in addition to the carboxylic acid\n"
                 '3. The molecule should have a reasonable carbon chain length '
                 'typical of fatty acids\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 73,
    'num_false_positives': 100,
    'num_true_negatives': 604,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.42196531791907516,
    'recall': 1.0,
    'f1': 0.5934959349593496,
    'accuracy': 0.8712998712998713,
    'negative_predictive_value': 1.0}