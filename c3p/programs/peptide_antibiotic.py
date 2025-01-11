"""
Classifies: CHEBI:25903 peptide antibiotic
"""
"""
Classifies: peptide antibiotic
"""
from rdkit import Chem

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    A peptide antibiotic is a peptide composed of multiple amino acid residues linked via peptide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a peptide antibiotic, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define peptide bond pattern (amide bond between amino acids)
    peptide_bond_pattern = Chem.MolFromSmarts("NC(=O)")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bonds)
    
    if num_peptide_bonds < 4:
        return False, f"Only {num_peptide_bonds} peptide bonds found, need at least 4"
    
    # Optionally, check for cyclic peptide (some peptide antibiotics are cyclic)
    ring_info = mol.GetRingInfo()
    has_macrocycle = any(len(ring) > 8 for ring in ring_info.AtomRings())
    if has_macrocycle:
        cyclic = "cyclic peptide"
    else:
        cyclic = "linear peptide"
    
    return True, f"Molecule is a {cyclic} with {num_peptide_bonds} peptide bonds"
    
__metadata__ = {
    'chemical_class': {
        'name': 'peptide antibiotic',
        'definition': 'A chemically diverse class of peptides that exhibit antimicrobial properties.'
    }
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25903',
                          'name': 'peptide antibiotic',
                          'definition': 'A chemically diverse class of '
                                        'peptides that exhibit antimicrobial '
                                        'properties.',
                          'parents': ['CHEBI:16670'],
                          'xrefs': [   'MetaCyc:Peptide-Antibiotics',
                                       'PMID:33398258'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O=C(NC1=C(C(=O)CCNC(=O)CCC(NC(=O)C)C(=O)OCC(O)CO)C=CC=C1)C(C=C)(C)C',
                                     'name': 'Citrinamide B',
                                     'reason': 'Only 3 peptide bonds found, '
                                               'need at least 4'},
                                 {   'smiles': 'S(CC[C@H](NC(=O)[C@@H](NC(=O)[C@@H](N)CC=1C=2C(NC1)=CC=CC2)CC(O)=O)C(O)=O)C',
                                     'name': 'Trp-Asp-Met',
                                     'reason': 'Only 2 peptide bonds found, '
                                               'need at least 4'},
                                 {   'smiles': 'C(CC)(O)(C)C',
                                     'name': '2-methylbutan-2-ol',
                                     'reason': 'Only 0 peptide bonds found, '
                                               'need at least 4'},
                                 {   'smiles': 'O1[C@@H]([C@H](O)[C@@H](O)[C@@H](NC(=O)C)C1O)CO',
                                     'name': 'N-Acetyl-D-Gulosamine',
                                     'reason': 'Only 1 peptide bonds found, '
                                               'need at least 4'},
                                 {   'smiles': '[C@@]12([C@]3([C@](CC[C@@]1([C@@]4(C(C[C@H](CC4)OS([O-])(=O)=O)=CC2)C)[H])([C@](CC3)([C@H](C)CC[C@@H](C(C)C)O)[H])C)[H])[H]',
                                     'name': '(24S)-hydroxycholesterol '
                                             '3-sulfate(1-)',
                                     'reason': 'Only 0 peptide bonds found, '
                                               'need at least 4'},
                                 {   'smiles': '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(C)O',
                                     'name': '17,20-dihydroxypregn-4-en-3-one',
                                     'reason': 'Only 0 peptide bonds found, '
                                               'need at least 4'},
                                 {   'smiles': 'CCCCCCCCCCc1c(C)c(O)c(OC)c(OC)c1O',
                                     'name': '6-decylubiquinol',
                                     'reason': 'Only 0 peptide bonds found, '
                                               'need at least 4'},
                                 {   'smiles': 'COC1=CC(=C(C=C1Br)OC)S(=O)(=O)N2CCCCCC2',
                                     'name': '1-(4-bromo-2,5-dimethoxyphenyl)sulfonylazepane',
                                     'reason': 'Only 0 peptide bonds found, '
                                               'need at least 4'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C=CC(=C2)N(C)C)O[C@@H]1CN(C)C(=O)NC3=CC=C(C=C3)F)[C@@H](C)CO',
                                     'name': '1-[[(2S,3R)-8-(dimethylamino)-5-[(2S)-1-hydroxypropan-2-yl]-3-methyl-6-oxo-3,4-dihydro-2H-1,5-benzoxazocin-2-yl]methyl]-3-(4-fluorophenyl)-1-methylurea',
                                     'reason': 'Only 3 peptide bonds found, '
                                               'need at least 4'},
                                 {   'smiles': 'P(=O)(N=C(N(CC)CC)C)(OCC)F',
                                     'name': 'A-234 nerve agent',
                                     'reason': 'Only 0 peptide bonds found, '
                                               'need at least 4'}],
    'sample_false_negatives': [   {   'smiles': 'COC1=CC(=O)N(C1)C(=O)[C@H]1CSC(=N1)C(CO)NC(=O)C1=CSC(\\C=N\\O)=N1',
                                      'name': 'althiomycin',
                                      'reason': 'Only 3 peptide bonds found, '
                                                'need at least 4'}],
    'attempt': 0,
    'reasoning': 'To classify peptide antibiotics based on their SMILES '
                 "string, we'll focus on identifying peptide structures with "
                 'multiple peptide bonds. Peptide antibiotics are peptides '
                 'composed of amino acid residues linked via peptide bonds '
                 '(amide bonds). \n'
                 '\n'
                 'The key steps are:\n'
                 '- Parse the SMILES string using RDKit.\n'
                 '- Identify peptide bonds (amide linkages) in the molecule.\n'
                 '- Count the number of peptide bonds.\n'
                 '- Set a threshold for the minimum number of peptide bonds '
                 '(e.g., at least 4) to consider the molecule as a peptide '
                 'antibiotic.\n'
                 '- Return True if the molecule meets the criteria, along with '
                 'a reason.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 37,
    'num_false_positives': 100,
    'num_true_negatives': 1784,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.27007299270072993,
    'recall': 0.9736842105263158,
    'f1': 0.4228571428571429,
    'accuracy': 0.9474505723204995,
    'negative_predictive_value': 0.9994397759103641}