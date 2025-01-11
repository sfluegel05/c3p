"""
Classifies: CHEBI:87659 dodecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester is any fatty acid ester where the carboxylic acid component is lauric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the lauric acid moiety pattern (C12 chain with a carboxyl group)
    lauric_acid_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC(=O)[OX2]")
    if not mol.HasSubstructMatch(lauric_acid_pattern):
        return False, "No lauric acid moiety found"

    # Check for ester bond (-COO-)
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) == 0:
        return False, "No ester bond found"

    # Verify that the ester bond is connected to the lauric acid moiety
    lauric_acid_match = mol.GetSubstructMatch(lauric_acid_pattern)
    lauric_acid_atoms = set(lauric_acid_match)
    ester_atoms = set()
    for match in ester_matches:
        ester_atoms.update(match)
    
    if not lauric_acid_atoms.intersection(ester_atoms):
        return False, "Ester bond not connected to lauric acid moiety"

    return True, "Contains lauric acid moiety with an ester bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:87659',
                          'name': 'dodecanoate ester',
                          'definition': 'Any fatty acid ester in which the '
                                        'carboxylic acid component is lauric '
                                        'acid.',
                          'parents': ['CHEBI:35748'],
                          'xrefs': ['PMID:23383323'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C1=CC=C2C(=C1)C3=C(C4=C(C=CC=N4)C=C3)NS2(=O)=O',
                                     'name': '5H-quinolino[8,7-c][1,2]benzothiazine '
                                             '6,6-dioxide',
                                     'reason': 'No lauric acid moiety found'},
                                 {   'smiles': 'C[C@@H]1CN([C@@H](COC2=C(C=C(C=C2)NS(=O)(=O)C3=CC=CC=C3)C(=O)N(C[C@H]1OC)C)C)C',
                                     'name': 'N-[(4R,7R,8S)-8-methoxy-4,5,7,10-tetramethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]benzenesulfonamide',
                                     'reason': 'No lauric acid moiety found'},
                                 {   'smiles': 'O=C1O[C@](C(=O)OC)(CC2=CC(=C(O)C=C2)CC(=O)C)C(=C1O)C3=CC=C(OC)C=C3',
                                     'name': 'Versicolactone A',
                                     'reason': 'No lauric acid moiety found'},
                                 {   'smiles': 'O=C(N[C@@H](CO)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC=1NC=NC1)CCC(O)=O',
                                     'name': 'His-Glu-Ser',
                                     'reason': 'No lauric acid moiety found'},
                                 {   'smiles': 'O=C1N(O)CCOCCNC(=O)CCC(=O)N(O)CCOCCNC(CCC(N(CCOCCNC(CC1)=O)O)=O)=O',
                                     'name': 'Desferrioxamine Et3',
                                     'reason': 'No lauric acid moiety found'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)(OCCN)(O)=O',
                                     'name': 'PE(18:3(9Z,12Z,15Z)/18:1(11Z))',
                                     'reason': 'No lauric acid moiety found'},
                                 {   'smiles': 'S(C[C@H](N)C(=O)N[C@H]1[C@@H](O)C(N2C3=NC=NC(=C3N=C2)N(C)C)OC1CO)C',
                                     'name': 'Cystocin',
                                     'reason': 'No lauric acid moiety found'},
                                 {   'smiles': 'COc1cc(O)c(CC=C(C)C)cc1C(=O)[C@@H](O)Cc1ccc(O)cc1',
                                     'name': 'lespeflorin C3',
                                     'reason': 'No lauric acid moiety found'},
                                 {   'smiles': 'NC(=O)OP(O)(O)=O',
                                     'name': 'carbamoyl phosphate',
                                     'reason': 'No lauric acid moiety found'},
                                 {   'smiles': 'O=C1OC(=CC(=C1C)O)C[C@@H](O)CC/C=C/C',
                                     'name': '6-(2′S-hydroxy-5′E-ene-1′-heptyl)-4-hydroxy-3-methyl-2H-pyran-2-one',
                                     'reason': 'No lauric acid moiety found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To classify chemical entities as dodecanoate esters, we need '
                 'to identify molecules that contain a lauric acid (dodecanoic '
                 'acid) component esterified to another molecule. Lauric acid '
                 'has the formula C12H24O2, and its esterification typically '
                 'involves the carboxyl group (-COOH) reacting with an alcohol '
                 'to form an ester bond (-COO-).\n'
                 '\n'
                 'The function `is_dodecanoate_ester` will check for the '
                 'presence of the lauric acid moiety (C12 chain with a '
                 'carboxyl group) and an ester bond. We will use RDKit to '
                 'parse the SMILES string and search for the relevant '
                 'substructures.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 65,
    'num_false_positives': 100,
    'num_true_negatives': 1382,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.3939393939393939,
    'recall': 1.0,
    'f1': 0.5652173913043478,
    'accuracy': 0.9353587588881707,
    'negative_predictive_value': 1.0}