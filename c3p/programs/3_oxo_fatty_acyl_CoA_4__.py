"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    This class features a 3-oxo group, a fatty acyl chain, and a Coenzyme A moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the Coenzyme A moiety pattern
    coa_pattern = Chem.MolFromSmarts("OP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"

    # Check for the 3-oxo group pattern
    oxo_group_pattern = Chem.MolFromSmarts("CC(=O)")
    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "No 3-oxo group found"

    # Check for long fatty acyl chain
    acyl_chain_pattern = Chem.MolFromSmarts("[#6]-[#6](=O)-[#6]-[#6]")
    if not mol.HasSubstructMatch(acyl_chain_pattern):
        return False, "No long fatty acyl chain found"

    return True, "Valid 3-oxo-fatty acyl-CoA(4-) structure identified"

# Example usage:
smiles = "CCCCCCCCCCCCCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
result, reason = is_3_oxo_fatty_acyl_CoA_4__(smiles)
print(f"Result: {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:57347',
                          'name': '3-oxo-fatty acyl-CoA(4-)',
                          'definition': 'An acyl-CoA(4-) arising from '
                                        'deprotonation of the phosphate and '
                                        'diphosphate groups of any 3-oxo-fatty '
                                        'acyl-CoA.',
                          'parents': ['CHEBI:77636', 'CHEBI:90726'],
                          'xrefs': ['MetaCyc:3-KETOACYL-COA'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O=C1C2=C(OC(=C1)C)C3=C(OC)C=C(OC)C=C3C(=C2O)C4=C5OC(=CC(C5=C(O)C=6C4=CC(OC)=CC6OC)=O)C',
                                     'name': 'Isonigerone',
                                     'reason': 'No Coenzyme A moiety found'},
                                 {   'smiles': 'CCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OC[C@@H](O)CO)OC(=O)CCCCCCC\\C=C/CCCCCCCC',
                                     'name': '1-myristoyl-2-oleoyl-sn-glycero-3-phosphatidylglycerol(1-)',
                                     'reason': 'No Coenzyme A moiety found'},
                                 {   'smiles': 'C(=C\\C/C=C\\CCCC(NC)=O)\\C/C=C\\C/C=C\\CCCCC',
                                     'name': 'N-methyl arachidonoyl amine',
                                     'reason': 'No Coenzyme A moiety found'},
                                 {   'smiles': 'OC1(O)[C@]23N(CC1)C(N(O)[C@H]([C@@]3(N=C(N2)N)[H])CO)=N',
                                     'name': 'Decarbamoylneosaxitoxin',
                                     'reason': 'No Coenzyme A moiety found'},
                                 {   'smiles': 'S([C@H]1N(C(=O)/C(=C\\C2=CC=C(OCC=C(C)C)C=C2)/N(C1=O)C)C)C',
                                     'name': 'Fusaperazine F',
                                     'reason': 'No Coenzyme A moiety found'},
                                 {   'smiles': 'Oc1ccccc1I',
                                     'name': '2-iodophenol',
                                     'reason': 'No Coenzyme A moiety found'},
                                 {   'smiles': 'O=C1N(CC(=O)N[C@H](C(=O)O[C@H]([C@@H](C(NC(C=C1)=C)=O)C)C(CCCCCCCCCCCCCC)C)C(O)C(=O)N)C',
                                     'name': 'Rakicidin H',
                                     'reason': 'No Coenzyme A moiety found'},
                                 {   'smiles': 'O(C1=C(OC)C=C(C2=C(O)C(OC)=C(C3=CC=CC=C3)C=C2OC)C=C1)C/C=C(/CO)\\C',
                                     'name': 'Prenylterphenyllin F',
                                     'reason': 'No Coenzyme A moiety found'},
                                 {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)[C@@H](O)[C@@H]1OC[C@H]3O[C@@H](OC[C@H]4O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]4O)[C@H](O)[C@@H](O)[C@@H]3O)CO[C@@H]5O[C@@H]([C@@H](O)[C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)[C@H]5O)CO[C@@H]7O[C@@H]([C@@H](O)[C@H](O)[C@H]7O)CO',
                                     'name': '(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5R,6R)-6-[[(2R,3R,4S,5R,6R)-3,5-Dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-2-yl]oxymethyl]-3,5-dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-3,4,5-trihydroxyoxan-2-yl]oxymethyl]oxane-2,3,4,5-tetrol',
                                     'reason': 'No Coenzyme A moiety found'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)O)[C@@H]4O[C@@H]([C@@H](O)[C@H](O[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O)CO)[C@H](O)[C@@H]6NC(=O)C)CO)CO)[C@@H]4O)CO[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO)[C@H](O)[C@H]9NC(=O)C)CO)CO',
                                     'name': 'Gal2GlcNAc2Man3GlcNAcFucGlcNAc',
                                     'reason': 'No Coenzyme A moiety found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To classify chemical entities of the class "3-oxo-fatty '
                 'acyl-CoA(4-)," we need to identify certain structural '
                 'features and patterns in the SMILES string, which are common '
                 'in all examples provided. The common features include:\n'
                 '\n'
                 '1. **Fatty acyl chain:** The SMILES should contain a '
                 'sequence of carbon atoms often interrupted by double bonds '
                 '(C=C), represented by segments like `\\C=C/` or `/C=C\\`. '
                 'This is indicative of unsaturated fatty acid chains.\n'
                 '   \n'
                 '2. **3-oxo group attachment:** The structure includes a '
                 'ketone group `(=O)` attached to an adjacent carbon, as '
                 'signified by `CC(=O)`.\n'
                 '\n'
                 '3. **Coenzyme A (CoA) moiety:** The presence of the CoA '
                 'structure, including the phosphate groups in the form '
                 '`OP([O-])(=O)` and the characteristic ring structure linked '
                 'to the phosphorus and nitrogen atoms. \n'
                 '\n'
                 'Given these observations, we should develop a function that '
                 'checks for these substructures in a given SMILES string and '
                 'confirms whether it corresponds to a "3-oxo-fatty '
                 'acyl-CoA(4-)."\n'
                 '\n'
                 'Below is the Python code to implement this classification:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 53,
    'num_false_positives': 83,
    'num_true_negatives': 142164,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.3897058823529412,
    'recall': 1.0,
    'f1': 0.5608465608465609,
    'accuracy': 0.9994167252283908,
    'negative_predictive_value': 1.0}