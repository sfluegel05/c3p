"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
from rdkit import Chem

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid.
    An olefinic fatty acid is characterized by having a carboxylic acid group
    and at least one carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an olefinic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Detect carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Detect carbon-carbon double bond (C=C)
    cc_double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(cc_double_bond_pattern):
        return False, "No carbon-carbon double bond found"

    return True, "Contains both a carboxylic acid group and at least one C=C double bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:53339',
                          'name': 'olefinic fatty acid',
                          'definition': 'Any fatty acid containing at least '
                                        'one C=C double bond.',
                          'parents': ['CHEBI:27208', 'CHEBI:78840'],
                          'xrefs': ['PMID:832335'],
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
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)CC4=CC=CC(=C4)C(=O)O)C5=CC=CC=C5N2C)[C@@H](C)CO',
                                     'name': 'LSM-9341',
                                     'reason': 'No carbon-carbon double bond '
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
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'OC(CCC1(CC)C(=O)NC(=O)NC1=O)(C)C',
                                      'name': "3'-Hydroxyamobarbital",
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'O=C1N(C(=O)C=2N=CN=NC2N1C)C',
                                      'name': 'Fervenulin',
                                      'reason': 'No carboxylic acid group '
                                                'found'}],
    'attempt': 0,
    'reasoning': 'To classify olefinic fatty acids, we need to design a '
                 'program that can identify the characteristic features of '
                 'molecules belonging to this class based on their SMILES '
                 'strings. \n'
                 '\n'
                 'Steps to classify olefinic fatty acids:\n'
                 '1. **Input Parsing**: Convert the given SMILES string into '
                 'an RDKit molecule object.\n'
                 '2. **Check for Validity**: Validate the molecule obtained '
                 'from the SMILES. If the SMILES is invalid, return `False` '
                 'along with a reason.\n'
                 '3. **Detect Carboxylic Acid Group**: An olefinic fatty acid '
                 'must contain a carboxylic acid group (`-COOH`). This can be '
                 'searched using a SMARTS pattern.\n'
                 '4. **Detect C=C Double Bonds**: Check the presence of at '
                 'least one C=C double bond. We use a SMARTS pattern to '
                 'identify this substructure.\n'
                 '5. **Return Classification**: If both features are present, '
                 'then the molecule is classified as an olefinic fatty acid.\n'
                 '\n'
                 'The following code implements this logic using RDKit:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 316,
    'num_false_positives': 100,
    'num_true_negatives': 576,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.7596153846153846,
    'recall': 0.9937106918238994,
    'f1': 0.8610354223433242,
    'accuracy': 0.89738430583501,
    'negative_predictive_value': 0.9965397923875432}