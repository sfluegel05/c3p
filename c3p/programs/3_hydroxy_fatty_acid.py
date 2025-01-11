"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: 3-hydroxy fatty acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
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
    
    # Look for beta-hydroxy pattern (OH on carbon 3 positions away from COOH)
    # Pattern: HO-C-C-C(=O)OH
    beta_hydroxy_pattern = Chem.MolFromSmarts("[OX2H1]-[CX4]-[CX4]-[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(beta_hydroxy_pattern):
        return False, "No hydroxy group at beta position (3-position)"

    # Count carbons to ensure it's a fatty acid (minimum 5 carbons for fatty acids)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 5:
        return False, "Carbon chain too short to be a fatty acid"
    
    # Count rotatable bonds to verify chain flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
        return False, "Too rigid to be a fatty acid"
    
    # Additional checks to exclude non-fatty acid structures
    # Check if molecule has reasonable number of oxygens (at least 3: COOH + OH)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 3:
        return False, "Insufficient oxygen atoms for 3-hydroxy fatty acid"
    
    # Verify the carbon chain is primarily linear (not too branched)
    # Count carbons with more than 2 carbon neighbors
    highly_branched_carbons = sum(1 for atom in mol.GetAtoms() 
                                if atom.GetAtomicNum() == 6 and 
                                len([n for n in atom.GetNeighbors() 
                                    if n.GetAtomicNum() == 6]) > 2)
    if highly_branched_carbons > carbon_count / 4:  # Allow some branching but not too much
        return False, "Carbon chain too branched for typical fatty acid"

    return True, "Contains carboxylic acid group and hydroxy group at 3-position with appropriate fatty acid chain"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59845',
                          'name': '3-hydroxy fatty acid',
                          'definition': 'Any fatty acid with a hydroxy '
                                        'functional group in the beta- or '
                                        '3-position. beta-Hydroxy fatty acids '
                                        'accumulate during cardiac hypoxia, '
                                        'and can also be used as chemical '
                                        'markers of bacterial endotoxins.',
                          'parents': ['CHEBI:24654', 'CHEBI:35969'],
                          'xrefs': [   'PMID:1226425',
                                       'PMID:17392575',
                                       'PMID:3094448',
                                       'PMID:7406063'],
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
                                 {   'smiles': 'S(O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]([C@@H](O)[C@H](O)CO[C@]2(O[C@H]([C@H](NC(=O)C)[C@@H](O)C2)[C@H](O)[C@H](O)CO)C(O)=O)[C@@H](NC(=O)C)CO)O[C@@H]1CO)(O)(=O)=O',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4R,5S)-5-acetamido-4-[(2R,3R,4R,5R,6R)-3,4-dihydroxy-6-(hydroxymethyl)-5-sulfooxyoxan-2-yl]oxy-2,3,6-trihydroxyhexoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No hydroxy group at beta '
                                               'position (3-position)'},
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
                                               'found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify 3-hydroxy fatty '
                 "acids. Based on the examples and definition, I'll identify "
                 'key characteristics:\n'
                 '\n'
                 '1. Must have a carboxylic acid group (-COOH)\n'
                 '2. Must have a hydroxy (-OH) group at the 3-position (beta '
                 'position)\n'
                 '3. Should have a carbon chain (fatty acid backbone)\n'
                 '4. May have additional hydroxy groups or other '
                 'modifications\n'
                 '5. May have specific stereochemistry (R/S) at the '
                 '3-position\n'
                 '\n'
                 "Here's the program:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 48,
    'num_false_positives': 100,
    'num_true_negatives': 5232,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.32432432432432434,
    'recall': 1.0,
    'f1': 0.489795918367347,
    'accuracy': 0.9814126394052045,
    'negative_predictive_value': 1.0}