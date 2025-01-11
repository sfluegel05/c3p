"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: B vitamin compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a B vitamin based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a B vitamin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for different B vitamins
    patterns = {
        # B1 (Thiamine) - thiazole ring connected to pyrimidine
        'B1': '[n+]1csc(CCO)c1C',
        
        # B2 (Riboflavin) - isoalloxazine ring system
        'B2': 'Cc1cc2nc3c(=O)[nX2]c(=O)nc-3n(CC(O)C(O)C(O)CO)c2cc1C',
        
        # B3 (Niacin/Nicotinic acid) - pyridine with carboxyl
        'B3': 'O=C(O)c1cccnc1',
        
        # B5 (Pantothenic acid) - beta-alanine with 2,4-dihydroxy-3,3-dimethylbutanamide
        'B5': 'CC(C)(CO)[CH](O)C(=O)NCCC(=O)[OH]',
        
        # B6 (Pyridoxine/Pyridoxal/Pyridoxamine) - pyridine with specific substitution
        'B6': 'Cc1ncc(CO)c(CO)c1O',
        
        # B7 (Biotin) - fused rings with sulfur
        'B7': 'O=C1NC2[CH]S[CH]C2N1',
        
        # B9 (Folate) - pterin + pABA + glutamate
        'B9': 'Nc1nc2ncc(CNc3ccc(CC)cc3)nc2c(=O)[nH]1',
        
        # B12 (Cobalamin) - corrin ring with cobalt
        'B12': '[Co]'
    }
    
    # Check each pattern
    for vitamin, pattern in patterns.items():
        substructure = Chem.MolFromSmarts(pattern)
        if substructure and mol.HasSubstructMatch(substructure):
            if vitamin == 'B12':
                # Additional check for corrin ring system for B12
                if any(atom.GetSymbol() == 'Co' for atom in mol.GetAtoms()):
                    return True, f"Matches vitamin {vitamin} (Cobalamin) structure"
            else:
                return True, f"Matches vitamin {vitamin} structure"
    
    # Additional checks for derivatives and variations
    
    # Check for phosphate groups (common in active forms)
    phosphate = Chem.MolFromSmarts('P(=O)(O)(O)O')
    
    # Check for nucleotide-like structures (as in FAD)
    nucleotide = Chem.MolFromSmarts('c1ncnc2[nH]cnc12')
    
    if mol.HasSubstructMatch(phosphate) and mol.HasSubstructMatch(nucleotide):
        return True, "Matches B vitamin derivative (likely FAD or FMN)"
    
    # Count number of rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    
    # Additional characteristics that might indicate a B vitamin
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if num_rings >= 2 and nitrogen_count >= 2 and oxygen_count >= 2:
        # This could be a modified or derivative form
        return True, "Matches structural characteristics of B vitamin derivative"
        
    return False, "Does not match any B vitamin structural patterns"


__metadata__ = {
    'chemical_class': {
        'name': 'B vitamin',
        'definition': 'Any member of the group of eight water-soluble vitamins '
                     'originally thought to be a single compound (vitamin B) that '
                     'play important roles in cell metabolism.',
        'subclasses': ['B1', 'B2', 'B3', 'B5', 'B6', 'B7', 'B9', 'B12']
    }
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:75769',
                          'name': 'B vitamin',
                          'definition': 'Any member of the group of eight '
                                        'water-soluble vitamins originally '
                                        'thought to be a single compound '
                                        '(vitamin B) that play important roles '
                                        'in cell metabolism. The group '
                                        'comprises of vitamin B1, B2, B3, B5, '
                                        'B6, B7, B9, and B12 (Around 20 other '
                                        'compounds were once thought to be B '
                                        'vitamins but are no longer classified '
                                        'as such).',
                          'parents': ['CHEBI:35352', 'CHEBI:36963'],
                          'xrefs': [   'MetaCyc:B-vitamins',
                                       'PMID:22743781',
                                       'PMID:23093174',
                                       'PMID:23238962',
                                       'PMID:23449527',
                                       'PMID:23462586',
                                       'PMID:23690582',
                                       'Wikipedia:B_vitamin'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'OC(=O)CCNC(O)=O',
                                     'name': 'N-carboxy-beta-alanine',
                                     'reason': 'Does not match any B vitamin '
                                               'structural patterns'},
                                 {   'smiles': 'O(C1=CC(=C(O)C=C1)C2=C(O)C=CC(=C2)OC)C',
                                     'name': '5,5′-dimethoxybiphenyl-2,2′-diol',
                                     'reason': 'Does not match any B vitamin '
                                               'structural patterns'},
                                 {   'smiles': 'COC1=CC=C(CCNC[C@H](O)C2=CC=C(O)C=C2)C=C1OC',
                                     'name': 'denopamine',
                                     'reason': 'Does not match any B vitamin '
                                               'structural patterns'},
                                 {   'smiles': '[Zn++].[S-]C(=S)NCCNC([S-])=S',
                                     'name': 'zineb',
                                     'reason': 'Does not match any B vitamin '
                                               'structural patterns'},
                                 {   'smiles': 'O=C1OC(O)C2=C1C[C@](O)([C@H]3CC(C[C@H]3[C@@H]2O)(C)C)C',
                                     'name': 'Lactarolide A',
                                     'reason': 'Does not match any B vitamin '
                                               'structural patterns'},
                                 {   'smiles': 'CC(=O)O[C@H]1CC[C@]23C[C@@H]1OO[C@@]2(C)C(=O)CCC3(C)C',
                                     'name': 'Talaperoxide B',
                                     'reason': 'Does not match any B vitamin '
                                               'structural patterns'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@H](OC(O)[C@@H]1NC(=O)C)CO[C@@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2O)CO)[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO',
                                     'name': 'N-[(3R,4R,5S,6R)-2,5-Dihydroxy-4-[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-3-yl]acetamide',
                                     'reason': 'Does not match any B vitamin '
                                               'structural patterns'},
                                 {   'smiles': 'COc1ccc(\\C=C/C2CCC=CC2c2ccc(OC)c(OC)c2)cc1OC',
                                     'name': '3-(3,4-dimethoxyphenyl)-4-[(Z)-3,4-dimethoxystyryl]cyclohex-1-ene',
                                     'reason': 'Does not match any B vitamin '
                                               'structural patterns'},
                                 {   'smiles': 'ClCC(P(=O)(O)O)=O',
                                     'name': 'Fosfonochlorin',
                                     'reason': 'Does not match any B vitamin '
                                               'structural patterns'},
                                 {   'smiles': '[H][C@]12C[C@@]3([H])O[C@@]4([H])CC=CCO[C@]4([H])[C@@H](O)[C@]3([H])O[C@]1([H])C=C[C@H](O)[C@@H](CO)O2',
                                     'name': '1,6:5,9:8,12:11,16-tetraanhydro-2,3,4,10,13,14-hexadeoxy-D-glycero-D-allo-D-gulo-heptadeca-2,13-dienitol',
                                     'reason': 'Does not match any B vitamin '
                                               'structural patterns'}],
    'sample_false_negatives': [   {   'smiles': '[H]C(=O)c1c(CO)cnc(C)c1O',
                                      'name': 'pyridoxal',
                                      'reason': 'Does not match any B vitamin '
                                                'structural patterns'},
                                  {   'smiles': 'Cc1ncc(CO)c(C[NH3+])c1O',
                                      'name': 'pyridoxaminium(1+)',
                                      'reason': 'Does not match any B vitamin '
                                                'structural patterns'},
                                  {   'smiles': '[H]C(=O)Cc1sc[n+](Cc2cnc(C)nc2N)c1C',
                                      'name': 'thiamine(1+) aldehyde',
                                      'reason': 'Does not match any B vitamin '
                                                'structural patterns'},
                                  {   'smiles': 'Cc1ncc(COP(O)(O)=O)c(CN)c1O',
                                      'name': "pyridoxamine 5'-phosphate",
                                      'reason': 'Does not match any B vitamin '
                                                'structural patterns'},
                                  {   'smiles': 'C1(O)=C(C)N=CC(CO)=C1C([H])=O.Cl',
                                      'name': 'pyridoxal hydrochloride',
                                      'reason': 'Does not match any B vitamin '
                                                'structural patterns'},
                                  {   'smiles': 'CC(C)(CO)[C@@H](O)C(=O)NCCC([O-])=O',
                                      'name': '(R)-pantothenate',
                                      'reason': 'Does not match any B vitamin '
                                                'structural patterns'},
                                  {   'smiles': 'C1(O)=C(C)[NH+]=CC(CO)=C1C[NH3+]',
                                      'name': 'pyridoxamine(2+)',
                                      'reason': 'Does not match any B vitamin '
                                                'structural patterns'},
                                  {   'smiles': '[H]C(=O)c1c(COP(O)(O)=O)cnc(C)c1O',
                                      'name': "pyridoxal 5'-phosphate",
                                      'reason': 'Does not match any B vitamin '
                                                'structural patterns'},
                                  {   'smiles': 'C1(O)=C(C)N=CC(CO)=C1CN.Cl',
                                      'name': 'pyridoxamine hydrochloride',
                                      'reason': 'Does not match any B vitamin '
                                                'structural patterns'},
                                  {   'smiles': 'Cc1ncc(CO)c(CN)c1O',
                                      'name': 'pyridoxamine',
                                      'reason': 'Does not match any B vitamin '
                                                'structural patterns'}],
    'attempt': 0,
    'reasoning': "I'll create a program to classify B vitamins. The key is to "
                 'identify the core structures of the different B vitamins '
                 '(B1/thiamine, B2/riboflavin, B3/niacin, B5/pantothenic acid, '
                 'B6/pyridoxine, B7/biotin, B9/folate, and B12/cobalamin).',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 72,
    'num_false_positives': 100,
    'num_true_negatives': 162,
    'num_false_negatives': 13,
    'num_negatives': None,
    'precision': 0.4186046511627907,
    'recall': 0.8470588235294118,
    'f1': 0.5603112840466926,
    'accuracy': 0.6743515850144092,
    'negative_predictive_value': 0.9257142857142857}