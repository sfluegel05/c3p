"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    Fatty acyl-CoA is characterized by having Coenzyme A moiety which includes a pattern with pantetheine unit and a nucleotide segment,
    as well as a fatty acyl chain typically attached via a thiol ester.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for key substructures indicating Coenzyme A
    # Search for common Coenzyme A substructure patterns
    coenzymeA_pattern = Chem.MolFromSmarts("NCCSC(=O)C")
    nucleotide_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)")

    # Check for presence of Coenzyme A segments
    if not mol.HasSubstructMatch(coenzymeA_pattern):
        return False, "Coenzyme A moiety pattern not found"
    
    if not mol.HasSubstructMatch(nucleotide_pattern):
        return False, "Nucleotide pattern not found"

    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate pattern not found"

    # Check for fatty acyl part via thiol ester and long carbon chains
    thiol_ester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thiol_ester_pattern):
        return False, "Thiol ester pattern not found"

    # Check for carbon chain length to verify fatty acyl part
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, f"Too few carbon atoms for a typical fatty acid chain"

    return True, "Contains structural features of a fatty acyl-CoA"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37554',
                          'name': 'fatty acyl-CoA',
                          'definition': 'An acyl-CoA that results from the '
                                        'formal condensation of the thiol '
                                        'group of coenzyme A with the carboxy '
                                        'group of any fatty acid.',
                          'parents': ['CHEBI:17984', 'CHEBI:61697'],
                          'xrefs': [   'PMID:11524729',
                                       'PMID:20442897',
                                       'PMID:2079609'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'Coenzyme A moiety pattern not '
                                               'found'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'Coenzyme A moiety pattern not '
                                               'found'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'Coenzyme A moiety pattern not '
                                               'found'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'Coenzyme A moiety pattern not '
                                               'found'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'Coenzyme A moiety pattern not '
                                               'found'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'Coenzyme A moiety pattern not '
                                               'found'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'Coenzyme A moiety pattern not '
                                               'found'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'Coenzyme A moiety pattern not '
                                               'found'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'Coenzyme A moiety pattern not '
                                               'found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'Coenzyme A moiety pattern not '
                                               'found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To classify a chemical as a fatty acyl-CoA, we need to '
                 'identify key structural features in its SMILES '
                 'representation. Fatty acyl-CoA molecules typically have the '
                 'following characteristics:\n'
                 '\n'
                 '1. **Pantoate Moiety**: Look for the bulky pantetheine unit '
                 'typically represented by segments like `NCCSC(=O)C`, '
                 'indicative of the presence of a "Coenzyme A" portion.\n'
                 '2. **Nucleotide Unit**: CoA contains a nucleotide-like '
                 'structure such as `n1cnc2c(N)ncnc12` toward the end.\n'
                 '3. **Phosphate Groups**: CoA often has phosphate groups, '
                 'represented by `OP(O)(=O)` patterns.\n'
                 '4. **Fatty Acyl Chain**: The fatty acyl group represented by '
                 'long carbon chains often terminated by a thiol ester '
                 '`C(=O)S`, which attaches to the CoA.\n'
                 '\n'
                 'Given these qualities, we can develop a function to detect '
                 'fatty acyl-CoA molecules by searching for this unique '
                 'combination of substructures within the SMILES string using '
                 'RDKIT.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 417,
    'num_false_positives': 100,
    'num_true_negatives': 14208,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.8065764023210832,
    'recall': 1.0,
    'f1': 0.892933618843683,
    'accuracy': 0.9932088285229203,
    'negative_predictive_value': 1.0}