"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: CHEBI:35341 3-oxo steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid has a steroid core structure with a ketone group at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid core pattern (four fused rings)
    steroid_core = Chem.MolFromSmarts("[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Look for ketone at position 3 
    # The pattern looks for a ketone group connected to the A ring of the steroid
    oxo_pattern = Chem.MolFromSmarts("[#6]1[#6][#6](=O)[#6][#6][#6]1")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No ketone group at position 3"

    # Additional checks to confirm steroid nature:
    
    # Count carbons (steroids typically have 17+ carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 17:
        return False, "Too few carbons for a steroid structure"

    # Check for proper ring connectivity
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings"

    # Look for common modifications that would invalidate the structure
    # (none identified that would specifically invalidate a 3-oxo steroid)

    return True, "Contains steroid core with ketone group at position 3"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35341',
        'name': '3-oxo steroid',
        'definition': 'Any oxo steroid where an oxo substituent is located at position 3.',
        'parents': ['CHEBI:35350']
    }
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:47788',
                          'name': '3-oxo steroid',
                          'definition': 'Any oxo steroid where an oxo '
                                        'substituent is located at position 3.',
                          'parents': ['CHEBI:35789', 'CHEBI:3992'],
                          'xrefs': [   'KEGG:C01876',
                                       'MetaCyc:3-Oxosteroids',
                                       'PMID:9811880'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(CC[C@]2(CC(C1)=O)[H])[H])(CC[C@@H]4O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)C(=O)O)O)O)O)[H])C)[H])C',
                                      'name': '5alpha-dihydrotestosterone '
                                              '17-O-(beta-D-glucuronide)',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC(=O)[C@@]4([H])CC(=O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)[C@@H](O)[C@H](O)[C@@H](C)C(C)C',
                                      'name': '3-dehydroteasterone',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])OC(C)=O',
                                      'name': 'testosterone acetate',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(=O)COC(=O)CCCCC',
                                      'name': 'hydrocortisone caproate',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@]12C[C@@H](O)[C@@]3([H])[C@]4([H])CC[C@]([H])([C@H](C)CCCC(C)C)[C@@]4(C)[C@@H](O)C[C@]3([H])[C@@]1(C)CCC(=O)C2',
                                      'name': '7alpha,12alpha-dihydroxy-5beta-cholestan-3-one',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@@]12CCC3=CC(=O)C=C[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])C[C@@H](C(=O)OC)[C@]1(O)C(=O)CO',
                                      'name': 'methyl '
                                              'prednisolone-16alpha-carboxylate',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3([C@@H](CC2=CC(C1)=O)O)[H])(CC[C@@]4([C@@H](CCC(=O)O)C)[H])[H])C)[H])C',
                                      'name': '7alpha-hydroxy-3-oxochol-4-en-24-oic '
                                              'acid',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'C[C@H](C[C@@H](O)C=C(C)C)[C@H]1CC[C@@]2(C)[C@@H]3[C@@H](O)C[C@@H]4[C@]5(C[C@@]35CC[C@]12C)CCC(=O)C4(C)C',
                                      'name': 'combretanone D',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': 'C[C@H]([C@H](O)[C@@H]1OC(=O)[C@@H](C)[C@@H]1C)[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@]4(C)[C@H]3CC[C@]12C',
                                      'name': 'paraminabeolide E',
                                      'reason': 'No steroid core structure '
                                                'found'},
                                  {   'smiles': '[H][C@@]12C[C@H](C)[C@](O)(C(=O)CO)[C@@]1(C)C[C@H](O)[C@@]1(Cl)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C',
                                      'name': 'beclomethasone',
                                      'reason': 'No steroid core structure '
                                                'found'}],
    'attempt': 0,
    'reasoning': "I'll create a program to classify 3-oxo steroids. Let's "
                 'break this down:\n'
                 '\n'
                 '1. A 3-oxo steroid has:\n'
                 '   - A steroid core (4 fused rings)\n'
                 '   - An oxo (=O) group at position 3\n'
                 '   - Various other substituents may be present\n'
                 '\n'
                 '2. Key characteristics:\n'
                 '   - Must have the basic steroid skeleton\n'
                 '   - Must have a ketone group at position 3\n'
                 '   - Can have various modifications but core structure must '
                 'remain\n'
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