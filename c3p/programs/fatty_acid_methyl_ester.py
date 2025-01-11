"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: fatty acid methyl ester
A fatty acid ester that is the carboxylic ester obtained by the formal condensation 
of a fatty acid with methanol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for methyl ester group pattern
    methyl_ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CH3X4]")
    methyl_ester_matches = mol.GetSubstructMatches(methyl_ester_pattern)
    
    if not methyl_ester_matches:
        return False, "No methyl ester group found"
    
    # Count number of methyl ester groups
    num_methyl_esters = len(methyl_ester_matches)
    
    # Special case: Allow dimethyl esters (like dimethyl sebacate)
    if num_methyl_esters > 2:
        return False, f"Too many methyl ester groups ({num_methyl_esters})"
    
    # Get atoms in methyl ester group
    ester_atoms = set()
    for match in methyl_ester_matches:
        ester_atoms.update(match)
    
    # Check for carbon chain attached to ester group
    carbon_chain = False
    for match in methyl_ester_matches:
        carbonyl_carbon = match[0]  # First atom in pattern is carbonyl carbon
        for neighbor in mol.GetAtomWithIdx(carbonyl_carbon).GetNeighbors():
            if neighbor.GetIdx() not in ester_atoms and neighbor.GetAtomicNum() == 6:
                # Found carbon attached to ester that's not part of ester group
                carbon_chain = True
                break
    
    if not carbon_chain:
        return False, "No carbon chain attached to ester group"
    
    # Count carbons (excluding methyl ester carbons)
    non_ester_carbons = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetIdx() not in ester_atoms:
            non_ester_carbons += 1
            
    if non_ester_carbons < 2:
        return False, "Carbon chain too short for fatty acid"
        
    # If we have exactly two methyl esters, check if it's a valid diester
    if num_methyl_esters == 2:
        # Check for linear chain between esters
        chain_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CH3X4].[CX3](=[OX1])[OX2][CH3X4]")
        if not mol.HasSubstructMatch(chain_pattern):
            return False, "Invalid diester structure"
        return True, "Valid dimethyl ester of dicarboxylic acid"
    
    return True, "Contains methyl ester group with appropriate carbon chain"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:4986',
                          'name': 'fatty acid methyl ester',
                          'definition': 'A fatty acid ester that is the '
                                        'carboxylic ester obtained by the '
                                        'formal condensation of a fatty acid '
                                        'with methanol.',
                          'parents': ['CHEBI:25248', 'CHEBI:35748'],
                          'xrefs': [   'KEGG:C03395',
                                       'MetaCyc:Fatty-acid-methyl-esters',
                                       'Wikipedia:Fatty_acid_methyl_ester'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C1C=C(N=C2N1NN=N2)C3=CC=CC=C3',
                                     'name': '5-phenyl-1,7-dihydrotetrazolo[1,5-a]pyrimidine',
                                     'reason': 'No methyl ester group found'},
                                 {   'smiles': 'O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)[C@H](COC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\CCCC',
                                     'name': 'TG(14:1(9Z)/22:1(13Z)/18:1(11Z))',
                                     'reason': 'No methyl ester group found'},
                                 {   'smiles': 'CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H]1O',
                                     'name': 'beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-{alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->6)]-alpha-D-Manp-(1->6)}-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)-beta-D-GlcpNAc',
                                     'reason': 'No methyl ester group found'},
                                 {   'smiles': '[C@@]12([C@@]([C@]([C@@H](CC1)C)(CCC(CC)C)C)(CCC[C@@H]2C)[H])C',
                                     'name': 'clerodane',
                                     'reason': 'No methyl ester group found'},
                                 {   'smiles': 'S(O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]([C@@H](O)[C@H](O)CO[C@]2(O[C@H]([C@H](NC(=O)C)[C@@H](O)C2)[C@H](O)[C@H](O)CO)C(O)=O)[C@@H](NC(=O)C)CO)O[C@@H]1CO)(O)(=O)=O',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4R,5S)-5-acetamido-4-[(2R,3R,4R,5R,6R)-3,4-dihydroxy-6-(hydroxymethyl)-5-sulfooxyoxan-2-yl]oxy-2,3,6-trihydroxyhexoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No methyl ester group found'},
                                 {   'smiles': 'O=C(NC(CC(C)C)C(=O)N[C@@H](CCCN=C(N)N)C=O)C(NC(=O)CC)CC(C)C',
                                     'name': 'Leupeptin Pr-LL',
                                     'reason': 'No methyl ester group found'},
                                 {   'smiles': 'COC1=CC=C(C=C1)CCNC(=O)C(=C2C3=CC=CC=C3C(=N2)NC(=O)C4=CC=CS4)C#N',
                                     'name': 'N-[3-[1-cyano-2-[2-(4-methoxyphenyl)ethylamino]-2-oxoethylidene]-1-isoindolyl]-2-thiophenecarboxamide',
                                     'reason': 'No methyl ester group found'},
                                 {   'smiles': 'O=C1O[C@@H](CC=2C1=C(OC)C(O)=CC2)CCC[C@@H](O)C',
                                     'name': 'Penicimarin C',
                                     'reason': 'No methyl ester group found'},
                                 {   'smiles': 'C[C@@H]1[C@H](O)CCC2=CC[C@H](C[C@]12C)C(C)=C',
                                     'name': '1-deoxycapsidiol',
                                     'reason': 'No methyl ester group found'},
                                 {   'smiles': 'C1CN=C(Cc2ccccc2)N1',
                                     'name': 'tolazoline',
                                     'reason': 'No methyl ester group found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify fatty acid methyl '
                 'esters (FAMEs). Based on the definition and examples, here '
                 'are the key characteristics:\n'
                 '\n'
                 '1. Must have a methyl ester group (-C(=O)OC)\n'
                 '2. Should have a carbon chain (fatty acid part)\n'
                 '3. Should not have other methyl esters (except in specific '
                 'cases like dimethyl esters)\n'
                 '4. Can have various modifications (unsaturations, epoxides, '
                 'hydroxy groups, etc.)\n'
                 '\n'
                 "Here's the program:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 141,
    'num_false_positives': 100,
    'num_true_negatives': 3476,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.5850622406639004,
    'recall': 1.0,
    'f1': 0.7382198952879581,
    'accuracy': 0.9730965832660748,
    'negative_predictive_value': 1.0}