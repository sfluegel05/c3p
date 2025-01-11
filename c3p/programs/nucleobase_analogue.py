"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: nucleobase analogue
"""

from rdkit import Chem

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is a molecule that can substitute for a normal nucleobase in nucleic acids.
    They typically contain pyrimidine-like or purine-like ring systems.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    has_pyrimidine_ring = False
    has_purine_ring = False

    # Get the ring information
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()

    for ring in atom_rings:
        ring_size = len(ring)
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]

        # Check if the ring is aromatic
        if all(atom.GetIsAromatic() for atom in ring_atoms):
            # For pyrimidine-like rings
            if ring_size == 6:
                nitrogen_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)
                if nitrogen_count >= 2:
                    has_pyrimidine_ring = True
            # For purine-like rings (fused rings)
            elif ring_size == 5 or ring_size == 6:
                # Check for fused rings
                bonds = [mol.GetBondBetweenAtoms(ring[i], ring[(i+1)%ring_size]) for i in range(ring_size)]
                bond_rings = [ri.NumBondRings(bond.GetIdx()) for bond in bonds]
                if max(bond_rings) > 1:
                    # Ring is part of a fused system
                    fused_atoms = set(ring)
                    for other_ring in atom_rings:
                        if other_ring != ring and len(set(other_ring).intersection(fused_atoms)) > 0:
                            # Found a fused ring system
                            combined_ring = set(ring).union(other_ring)
                            combined_ring_atoms = [mol.GetAtomWithIdx(idx) for idx in combined_ring]
                            nitrogen_count = sum(1 for atom in combined_ring_atoms if atom.GetAtomicNum() == 7)
                            if len(combined_ring) >= 9 and nitrogen_count >= 4:
                                has_purine_ring = True
                                break
    if has_purine_ring:
        return True, "Contains purine-like fused ring system"
    elif has_pyrimidine_ring:
        return True, "Contains pyrimidine-like ring system"
    else:
        return False, "Does not contain purine or pyrimidine-like ring systems"
        

__metadata__ = {   'chemical_class': {   'id': None,
                              'name': 'nucleobase analogue',
                              'definition': 'A molecule that can substitute for a normal nucleobase in nucleic acids.',
                              'parents': []},
        'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8,
                      'max_attempts': 5,
                      'max_positive_instances': None,
                      'max_positive_to_test': None,
                      'max_negative_to_test': None,
                      'max_positive_in_prompt': 50,
                      'max_negative_in_prompt': 20,
                      'max_instances_in_prompt': 100,
                      'test_proportion': 0.1},
        'message': None,
        'attempt': 0,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:67142',
                          'name': 'nucleobase analogue',
                          'definition': 'A molecule that can substitute for a '
                                        'normal nucleobase in nucleic acids.',
                          'parents': ['CHEBI:33832'],
                          'xrefs': ['Wikipedia:Base_analog'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'CNC(=O)CC[C@H](N)C(O)=O',
                                     'name': 'N(5)-methyl-L-glutamine',
                                     'reason': 'Does not contain purine or '
                                               'pyrimidine-like ring systems'},
                                 {   'smiles': '[H][C@@]1(CC[C@]2(C)[C@]1([H])CC[C@]1([H])[C@@]3(C)CCCC(C)(C)[C@]3([H])CC[C@@]21C)C(=C)CCC=C(C)C',
                                     'name': 'dammara-20,24-diene',
                                     'reason': 'Does not contain purine or '
                                               'pyrimidine-like ring systems'},
                                 {   'smiles': 'O=C1C2=C(O)C3=C(O)C=CC=C3C(=C2C(C(=O)C)C(C1)(O)C)C=4C=5C(C(O)=C6C4C(C(=O)C)C(O)(C)CC6=O)=C(O)C=CC5',
                                     'name': 'A 39183A',
                                     'reason': 'Does not contain purine or '
                                               'pyrimidine-like ring systems'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@H]1OC)C)C)CC3=CC=C(C=C3)F',
                                     'name': 'N-[(4S,7R,8S)-5-[(4-fluorophenyl)methyl]-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'Does not contain purine or '
                                               'pyrimidine-like ring systems'},
                                 {   'smiles': 'C(CCN1CCCCC1)(CCN(C(C)C)C(C)=O)(C(N)=O)C2=C(Cl)C=CC=C2',
                                     'name': 'bidisomide',
                                     'reason': 'Does not contain purine or '
                                               'pyrimidine-like ring systems'},
                                 {   'smiles': 'OCCCCCCCCCCC#CC=C',
                                     'name': '13-Tetradece-11-yn-1-ol',
                                     'reason': 'Does not contain purine or '
                                               'pyrimidine-like ring systems'},
                                 {   'smiles': 'C[C@@H]1CCCCO[C@@H]([C@@H](CN(C(=O)C2=C(O1)C=CC(=C2)N(C)C)[C@@H](C)CO)C)CN(C)CC3=CC=C(C=C3)OC',
                                     'name': '(3R,9S,10R)-16-(dimethylamino)-12-[(2S)-1-hydroxypropan-2-yl]-9-[[(4-methoxyphenyl)methyl-methylamino]methyl]-3,10-dimethyl-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-13-one',
                                     'reason': 'Does not contain purine or '
                                               'pyrimidine-like ring systems'},
                                 {   'smiles': 'O(CCC12CC3CC(C1)CC(C2)C3)C(=O)CC4=CC=C(OCC(O)CNC(C)C)C=C4',
                                     'name': 'Adaprolol',
                                     'reason': 'Does not contain purine or '
                                               'pyrimidine-like ring systems'},
                                 {   'smiles': 'CC1=CC=CC=C1C2=CC=C(C=C2)C3[C@H]4CNC[C@@H]3N4C(=O)NC5=CC=C(C=C5)Cl',
                                     'name': '(1S,5R)-N-(4-chlorophenyl)-7-[4-(2-methylphenyl)phenyl]-3,6-diazabicyclo[3.1.1]heptane-6-carboxamide',
                                     'reason': 'Does not contain purine or '
                                               'pyrimidine-like ring systems'},
                                 {   'smiles': 'O([C@@H]1[C@@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)OC[C@@H](O)[C@H](O)[C@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO)[C@@H](NC(=O)C)CO)O[C@@H]([C@@H]1O)CO)[C@@H]4O[C@@H]([C@@H](O[C@@H]5O[C@H]([C@@H](O)[C@@H](O)[C@@H]5O)C)[C@H](O)[C@H]4NC(=O)C)CO',
                                     'name': 'N-[(2R,3R,4R,5S,6R)-5-[(2S,3R,4S,5S,6R)-4-[(2S,3R,4R,5S,6R)-3-Acetamido-4-hydroxy-6-(hydroxymethyl)-5-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2-[(2R,3S,4R,5S)-5-acetamido-2,3,6-trihydroxy-4-[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyhexoxy]-4-hydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'Does not contain purine or '
                                               'pyrimidine-like ring systems'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To determine if a molecule is a nucleobase analogue, we can '
                 'check if it contains pyrimidine-like or purine-like ring '
                 'systems:\n'
                 '\n'
                 '- **Pyrimidine-like rings**: aromatic six-membered rings '
                 'with at least two nitrogen atoms.\n'
                 '- **Purine-like rings**: fused aromatic ring systems '
                 'consisting of a five-membered and a six-membered ring with '
                 'at least four nitrogen atoms in total.\n'
                 '\n'
                 'Below is the code that implements this logic:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 37,
    'num_false_positives': 100,
    'num_true_negatives': 1597,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.27007299270072993,
    'recall': 1.0,
    'f1': 0.42528735632183906,
    'accuracy': 0.9423298731257209,
    'negative_predictive_value': 1.0}