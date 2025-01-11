"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: methyl sulfide
Definition: Any aliphatic sulfide in which at least one of the organyl groups attached to the sulfur is a methyl group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule contains a methyl sulfide group based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for S-CH3 pattern
    # [#16X2] = divalent sulfur
    # [CH3X4] = sp3 methyl group
    methylsulfide_pattern = Chem.MolFromSmarts("[#16X2][CH3X4]")
    matches = mol.GetSubstructMatches(methylsulfide_pattern)
    
    if not matches:
        return False, "No methylsulfide (S-CH3) group found"
    
    # Check that sulfur is not part of a higher oxidation state group
    # Look for S=O, S(=O)=O patterns that would indicate sulfoxide/sulfone
    sulfoxide_pattern = Chem.MolFromSmarts("[#16X3](=[OX1])")
    sulfone_pattern = Chem.MolFromSmarts("[#16X4](=[OX1])=[OX1]")
    
    for match in matches:
        s_atom_idx = match[0]  # Index of the sulfur atom
        s_atom = mol.GetAtomWithIdx(s_atom_idx)
        
        # Check if this sulfur is part of sulfoxide/sulfone
        if mol.HasSubstructMatch(sulfoxide_pattern) or mol.HasSubstructMatch(sulfone_pattern):
            continue
            
        # Verify sulfur has exactly two single bonds (sulfide)
        if s_atom.GetDegree() == 2 and s_atom.GetTotalValence() == 2:
            return True, "Contains methylsulfide (S-CH3) group"
            
    return False, "Sulfur is not in sulfide form (may be sulfoxide, sulfone, or other S species)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:86315',
                          'name': 'methyl sulfide',
                          'definition': 'Any aliphatic sulfide in which at '
                                        'least one of the organyl groups '
                                        'attached to the sulfur is a methyl '
                                        'group.',
                          'parents': ['CHEBI:16385'],
                          'xrefs': ['MetaCyc:Methyl-thioethers'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O1[C@H](CC2=C(C=3C=4C(C(O)=C5C3C[C@H](C)OC5)=C(OC)C=C(C4)OC)C6=CC(OC)=CC(=C6C(=C2C1)O)OC)C',
                                     'name': '(3S)-5-[(3S)-10-hydroxy-7,9-dimethoxy-3-methyl-3,4-dihydro-1H-benzo[g]isochromen-5-yl]-7,9-dimethoxy-3-methyl-3,4-dihydro-1H-benzo[g]isochromen-10-ol',
                                     'reason': 'No methylsulfide (S-CH3) group '
                                               'found'},
                                 {   'smiles': 'C[C@H]1CN([C@@H](COC2=C(C=CC(=C2)NC(=O)C)C(=O)N(C[C@@H]1OC)C)C)C(=O)C3CCC3',
                                     'name': 'N-[(5R,6S,9R)-8-[cyclobutyl(oxo)methyl]-5-methoxy-3,6,9-trimethyl-2-oxo-11-oxa-3,8-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'No methylsulfide (S-CH3) group '
                                               'found'},
                                 {   'smiles': '[C@@]12(CCCC[C@@]1(C=C[C@@H]([C@@H]2CC[C@H](C[C@H](CC(=O)[O-])O)O)C)[H])[H]',
                                     'name': '4a,5-dihydro-ML-236C carboxylate',
                                     'reason': 'No methylsulfide (S-CH3) group '
                                               'found'},
                                 {   'smiles': 'S(=O)(CC1OC(N2C3=NC=NC(N)=C3N=C2)C(O)C1O)C',
                                     'name': "(S)-5'-Deoxy-5'-(methylsulfinyl)adenosine",
                                     'reason': 'No methylsulfide (S-CH3) group '
                                               'found'},
                                 {   'smiles': 'OC=1C(=C(C=2C=3C(NC2)=CC=CC3)C(=O)C(=O)C1C=4C=5C(NC4)=CC=CC5)CC=C(C)C',
                                     'name': 'Ochrindole D',
                                     'reason': 'No methylsulfide (S-CH3) group '
                                               'found'},
                                 {   'smiles': 'O(C1C(O)C(OC1OC2=C(OC=3C(C2=O)=C(O)C=C(OC4OC(C(O)C(O)C4O)C)C3)C5=CC=C(O)C=C5)CO)C6OCC(O)(C6O)CO',
                                     'name': 'Kaempferol '
                                             '3-apiosyl-(1->2)-alpha-L-arabinofuranoside-7-rhamnoside',
                                     'reason': 'No methylsulfide (S-CH3) group '
                                               'found'},
                                 {   'smiles': '[O-]C(=O)OON=O',
                                     'name': 'nitrosoperoxycarbonate(1-)',
                                     'reason': 'No methylsulfide (S-CH3) group '
                                               'found'},
                                 {   'smiles': 'P(OCC(O)COC(=O)CCCCCCCCCCCCCCCCC)(O)(O)=O',
                                     'name': 'LysoPA(18:0/0:0)',
                                     'reason': 'No methylsulfide (S-CH3) group '
                                               'found'},
                                 {   'smiles': 'O=C1C=2C(=O)OC34C2OC(C(O)C3C=CC(C(CC1C)C)=O)C(C)C4',
                                     'name': 'Atrop-Abybetaomicin C',
                                     'reason': 'No methylsulfide (S-CH3) group '
                                               'found'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCCCCC)CO/C=C\\CCCCCCCCCCCCCCCC)(OC[C@@H](O)CO)(O)=O',
                                     'name': 'PG(P-18:0/19:0)',
                                     'reason': 'No methylsulfide (S-CH3) group '
                                               'found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify methyl sulfides. '
                 'The key characteristics are:\n'
                 '1. Must contain at least one S-CH3 group '
                 '(methylsulfanyl/methylthio group)\n'
                 '2. The sulfur must be connected via single bonds (sulfide)\n'
                 "3. No restrictions on what's on the other side of the "
                 'sulfur\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 69,
    'num_false_positives': 100,
    'num_true_negatives': 6848,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.40828402366863903,
    'recall': 1.0,
    'f1': 0.5798319327731092,
    'accuracy': 0.9857488955394043,
    'negative_predictive_value': 1.0}