"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: branched-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight (CoA itself is quite large)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 750:  # CoA MW is ~767
        return False, "Molecular weight too low for acyl-CoA"

    # Count key atoms
    num_P = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15])  # Phosphorus
    num_S = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16])  # Sulfur
    num_N = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7])   # Nitrogen
    
    # Basic CoA checks
    if num_P < 3:
        return False, "Missing phosphate groups characteristic of CoA"
    if num_S < 1:
        return False, "Missing sulfur atom required for thioester linkage"
    if num_N < 5:
        return False, "Insufficient nitrogen atoms for CoA structure"

    # Check for key CoA structural features
    coenzyme_a_patterns = [
        # Thioester linkage
        "[CX3](=O)[SX2]",
        # Adenine base
        "n1c(nc2c1ncnc2N)",
        # Pantetheine part
        "CC(C)(COP(=O)(O)OP(=O)(O)O)C"
    ]
    
    for pattern in coenzyme_a_patterns:
        if not mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, "Missing essential CoA structural feature"

    # Check for branching in the fatty acid portion
    branching_patterns = [
        # General branched carbon
        "[CH1,CH0](-[CH3])(-[CH2,CH3])-[!R]",
        # Terminal isopropyl group
        "[CH3]-[CH](-[CH3])-[CH2]",
        # Internal branching
        "[CH2]-[CH](-[CH3])-[CH2]"
    ]
    
    found_branching = False
    for pattern in branching_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_branching = True
            break
            
    if not found_branching:
        return False, "No branching pattern found in fatty acid portion"

    # Check for fatty acid chain length (at least 4 carbons)
    carbon_chain = "[CH2,CH1,CH0]-[CH2,CH1,CH0]-[CH2,CH1,CH0]-[CH2,CH1,CH0]"
    if not mol.HasSubstructMatch(Chem.MolFromSmarts(carbon_chain)):
        return False, "Carbon chain too short for fatty acid"

    # Look for common modifications
    modifications = []
    
    # Check for hydroxy groups
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CH1,CH0]([OH1])")):
        modifications.append("hydroxy-modified")
    
    # Check for unsaturation
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C=C")):
        modifications.append("unsaturated")
    
    # Check for oxo groups (besides thioester)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CH2]C(=O)[CH2]")):
        modifications.append("oxo-modified")

    base_reason = "Contains CoA moiety with thioester linkage and branched-chain fatty acid pattern"
    if modifications:
        return True, f"{base_reason} ({', '.join(modifications)})"
    return True, base_reason


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61912',
                          'name': 'branched-chain fatty acyl-CoA',
                          'definition': 'A fatty acyl-CoA that results from '
                                        'the formal condensation of the thiol '
                                        'group of coenzyme A with the carboxy '
                                        'group of any branched-chain fatty '
                                        'acid.',
                          'parents': ['CHEBI:37554'],
                          'xrefs': [   'PMID:16768463',
                                       'PMID:1761563',
                                       'PMID:9553048'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               "Error: 'Mol' object has no attribute 'HasSubstructMatches'\n"
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C1C=C(N=C2N1NN=N2)C3=CC=CC=C3',
                                     'name': '5-phenyl-1,7-dihydrotetrazolo[1,5-a]pyrimidine',
                                     'reason': 'Molecular weight too low for '
                                               'acyl-CoA'},
                                 {   'smiles': 'O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)[C@H](COC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\CCCC',
                                     'name': 'TG(14:1(9Z)/22:1(13Z)/18:1(11Z))',
                                     'reason': 'Missing phosphate groups '
                                               'characteristic of CoA'},
                                 {   'smiles': 'CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H]1O',
                                     'name': 'beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-{alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->6)]-alpha-D-Manp-(1->6)}-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)-beta-D-GlcpNAc',
                                     'reason': 'Missing phosphate groups '
                                               'characteristic of CoA'},
                                 {   'smiles': '[C@@]12([C@@]([C@]([C@@H](CC1)C)(CCC(CC)C)C)(CCC[C@@H]2C)[H])C',
                                     'name': 'clerodane',
                                     'reason': 'Molecular weight too low for '
                                               'acyl-CoA'},
                                 {   'smiles': 'S(O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]([C@@H](O)[C@H](O)CO[C@]2(O[C@H]([C@H](NC(=O)C)[C@@H](O)C2)[C@H](O)[C@H](O)CO)C(O)=O)[C@@H](NC(=O)C)CO)O[C@@H]1CO)(O)(=O)=O',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4R,5S)-5-acetamido-4-[(2R,3R,4R,5R,6R)-3,4-dihydroxy-6-(hydroxymethyl)-5-sulfooxyoxan-2-yl]oxy-2,3,6-trihydroxyhexoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'Missing phosphate groups '
                                               'characteristic of CoA'},
                                 {   'smiles': 'O=C(NC(CC(C)C)C(=O)N[C@@H](CCCN=C(N)N)C=O)C(NC(=O)CC)CC(C)C',
                                     'name': 'Leupeptin Pr-LL',
                                     'reason': 'Molecular weight too low for '
                                               'acyl-CoA'},
                                 {   'smiles': 'COC1=CC=C(C=C1)CCNC(=O)C(=C2C3=CC=CC=C3C(=N2)NC(=O)C4=CC=CS4)C#N',
                                     'name': 'N-[3-[1-cyano-2-[2-(4-methoxyphenyl)ethylamino]-2-oxoethylidene]-1-isoindolyl]-2-thiophenecarboxamide',
                                     'reason': 'Molecular weight too low for '
                                               'acyl-CoA'},
                                 {   'smiles': 'O=C1O[C@@H](CC=2C1=C(OC)C(O)=CC2)CCC[C@@H](O)C',
                                     'name': 'Penicimarin C',
                                     'reason': 'Molecular weight too low for '
                                               'acyl-CoA'},
                                 {   'smiles': 'CO[C@@H]1[C@@H]2[C@@H](C[C@@H]3[C@@](O2)([C@H]([C@@H]([C@H](O3)C(=O)OC)O)O)O)OC1N4C=NC5=C4N=CN=C5N',
                                     'name': 'LSM-4497',
                                     'reason': 'Molecular weight too low for '
                                               'acyl-CoA'},
                                 {   'smiles': 'C[C@@H]1[C@H](O)CCC2=CC[C@H](C[C@]12C)C(C)=C',
                                     'name': '1-deoxycapsidiol',
                                     'reason': 'Molecular weight too low for '
                                               'acyl-CoA'}],
    'sample_false_negatives': [],
    'attempt': 3,
    'reasoning': 'Let me analyze the previous code and suggest improvements:\n'
                 '\n'
                 '1. Error Analysis:\n'
                 "- The main error 'Mol' object has no attribute "
                 "'HasSubstructMatches' suggests we used the wrong method "
                 'name. It should be `HasSubstructMatch` (singular) not '
                 '`HasSubstructMatches` (plural).\n'
                 '- The code was trying to be too specific with SMARTS '
                 'patterns that may not match all variations.\n'
                 '- The branching patterns may have been too restrictive and '
                 'not covering all possible cases.\n'
                 '\n'
                 '2. Improvements needed:\n'
                 '- Use correct method names for RDKit\n'
                 '- Simplify the branching detection logic\n'
                 '- Make the CoA pattern detection more robust\n'
                 '- Add molecular weight and atom count checks\n'
                 '- Consider using simpler SMARTS patterns\n'
                 '- Add better handling of unsaturated bonds\n'
                 '\n'
                 "Here's my improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 45,
    'num_false_positives': 100,
    'num_true_negatives': 9485,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.3103448275862069,
    'recall': 1.0,
    'f1': 0.4736842105263158,
    'accuracy': 0.9896157840083074,
    'negative_predictive_value': 1.0}