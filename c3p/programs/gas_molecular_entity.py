"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: CHEBI:33262 gas molecular entity
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is likely to be gaseous at STP (0°C and 100 kPa)
    based on its molecular structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a gas at STP, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get basic molecular properties
    mol_wt = Descriptors.ExactMolWt(mol)
    num_atoms = mol.GetNumAtoms()
    num_heavy_atoms = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    
    # Check for noble gases (He, Ne, Ar, Kr, Xe, Rn)
    noble_gases = {'He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn'}
    if num_atoms == 1:
        symbol = mol.GetAtomWithIdx(0).GetSymbol()
        if symbol in noble_gases:
            return True, f"Noble gas ({symbol})"
        if symbol == 'C':
            return True, "Carbon atom"
            
    # Check molecular weight - most gases are light
    if mol_wt > 150:  # Most common gases are below this weight
        # Special case for perfluorinated compounds which can be gases at higher weights
        num_F = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[F]')))
        if num_F < 6:  # Allow higher weight only for heavily fluorinated compounds
            return False, f"Molecular weight too high ({mol_wt:.1f} Da)"
    
    # Most gases have relatively few heavy atoms
    if num_heavy_atoms > 8 and 'F' not in smiles:
        return False, f"Too many heavy atoms ({num_heavy_atoms})"
        
    # Common gas patterns
    gas_patterns = {
        'diatomic': '[!C!H]~[!C!H]',  # For N2, O2, F2, Cl2, etc.
        'carbon_oxide': '[C]=[O,N]',  # For CO, CO2
        'small_alkane': '[CX4]~[CX4]',  # For CH4, C2H6, etc.
        'small_alkene': '[CX3]=[CX3]',  # For C2H4, etc.
        'small_alkyne': '[CX2]#[CX2]',  # For C2H2, etc.
        'hydrogen_halide': '[F,Cl,Br,I][H]',  # For HF, HCl, HBr, HI
        'ammonia': '[NX3]([H])([H])[H]',  # NH3
        'ozone': '[O-][O+]=O'  # O3
    }
    
    for name, pattern in gas_patterns.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return True, f"Matches {name} pattern"
            
    # Check for small hydrocarbons (up to C4)
    if all(atom.GetSymbol() in ('C', 'H') for atom in mol.GetAtoms()):
        if num_heavy_atoms <= 4:
            return True, f"Small hydrocarbon with {num_heavy_atoms} carbons"
            
    # Check for simple molecules with few non-H atoms
    if num_heavy_atoms <= 3:
        return True, f"Simple molecule with {num_heavy_atoms} heavy atoms"
        
    # If we haven't returned True by now, it's probably not a gas
    return False, "Does not match known gas patterns and is too complex"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:33262',
        'name': 'gas molecular entity',
        'definition': 'Any main group molecular entity that is gaseous at standard temperature and pressure (STP; 0degreeC and 100 kPa).',
        'parents': []
    }
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:138675',
                          'name': 'gas molecular entity',
                          'definition': 'Any main group molecular entity that '
                                        'is gaseous at standard temperature '
                                        'and pressure (STP; 0degreeC and 100 '
                                        'kPa).',
                          'parents': ['CHEBI:33579'],
                          'xrefs': [   'Wikipedia:https://en.wikipedia.org/wiki/Gas'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O=C1O[C@H]2[C@H](O[C@H](C2)C[C@H](O)C)C=3C1=C(O)C(OC)=C(OC)C3',
                                     'name': '(12R)-12-hydroxymonocerin',
                                     'reason': 'Molecular weight too high '
                                               '(324.1 Da)'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCC/C=C\\CCCCCCCC)(OCCN)(O)=O',
                                     'name': 'PE(20:1(11Z)/16:0)',
                                     'reason': 'Molecular weight too high '
                                               '(745.6 Da)'},
                                 {   'smiles': 'O=C(C1=C(N)C2=C(OC(C)(C)C=C2)C=C1)C[C@@H](NC(=O)C)CO',
                                     'name': 'Fusarochromene',
                                     'reason': 'Molecular weight too high '
                                               '(318.2 Da)'},
                                 {   'smiles': 'COc1cc(\\C=C\\C(=O)c2c(O)cc(O)cc2O)ccc1O',
                                     'name': 'homoeriodictyol chalcone',
                                     'reason': 'Molecular weight too high '
                                               '(302.1 Da)'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@H](O[C@@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@@H]([C@H](O)[C@@H](NC(=O)C)CO)[C@H](O)CO)[C@H]1O)CO[C@H]3O[C@@H]([C@@H](O)[C@H](O)[C@@H]3O[C@@H]4O[C@@H]([C@@H](O[C@@H]5O[C@@H]([C@H](O)[C@H](O[C@]6(O[C@H]([C@H](NC(=O)C)[C@@H](O)C6)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]5O)CO)[C@H](O)[C@H]4NC(=O)C)CO)CO)[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@@H]7O[C@@H]8O[C@@H]([C@@H](O[C@@H]9O[C@@H]([C@H](O)[C@H](O)[C@H]9NC(=O)C)CO[C@]%10(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%10)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]8NC(=O)C)CO)CO',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[[(2R,3R,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4S,5S,6R)-2-[(2R,3R,4S,5S,6S)-2-[[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-5-[(2S,3R,4S,5S,6R)-4-[(2S,4S,5R,6R)-5-acetamido-2-carboxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxan-2-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5S)-5-acetamido-1,2,4,6-tetrahydroxyhexan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,5-dihydroxyoxan-4-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4-dihydroxyoxan-2-yl]methoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'Molecular weight too high '
                                               '(2265.8 Da)'},
                                 {   'smiles': 'O1C=2C(CC=C(C)C)=C(OC)C=C(O)C2C(=O)C=C1C3=CC=C(O)C=C3',
                                     'name': 'Artonin U',
                                     'reason': 'Molecular weight too high '
                                               '(352.1 Da)'},
                                 {   'smiles': 'O=C(CCCCCCCCC)CCCCCC(O)=O',
                                     'name': '7-Keto palmitic acid',
                                     'reason': 'Molecular weight too high '
                                               '(270.2 Da)'},
                                 {   'smiles': 'CNC(=N)NCCC[C@H]([NH3+])C([O-])=O',
                                     'name': 'N(omega)-methyl-L-arginine '
                                             'zwitterion',
                                     'reason': 'Molecular weight too high '
                                               '(188.1 Da)'},
                                 {   'smiles': '[C@]([C@@](/C=C/CCCCCCCCCCC)(O)[H])(NC(=O)C(CCCCCCCCCCCCCCCCCCCCCC)O)([H])CO',
                                     'name': 'N-(2-hydroxy-tetracosanoyl)-hexadecasphing-4-enine',
                                     'reason': 'Molecular weight too high '
                                               '(637.6 Da)'},
                                 {   'smiles': 'COc1cc2oc(-c3ccc(O)cc3O)c(CC=C(C)C)c(=O)c2c(O)c1\\C=C\\C(C)C',
                                     'name': 'artocarpin',
                                     'reason': 'Molecular weight too high '
                                               '(436.2 Da)'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify gas molecular '
                 'entities. Based on the examples, I notice that these are '
                 'generally:\n'
                 '\n'
                 '1. Noble gases (He, Ne, Ar, Kr, Xe, Rn)\n'
                 '2. Small molecules (typically <5 non-H atoms)\n'
                 '3. Simple hydrocarbons (alkanes, alkenes, alkynes)\n'
                 '4. Simple inorganic molecules (CO, CO2, O2, O3, NH3, HCl, '
                 'HI)\n'
                 '\n'
                 'The key characteristics appear to be:\n'
                 '- Low molecular weight (generally <100 Da)\n'
                 '- Simple structure (few atoms, simple bonds)\n'
                 '- Limited branching\n'
                 '- Often symmetric\n'
                 '- Can include isotopes\n'
                 '\n'
                 "Here's my implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 42,
    'num_false_positives': 100,
    'num_true_negatives': 4179,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.29577464788732394,
    'recall': 1.0,
    'f1': 0.45652173913043476,
    'accuracy': 0.9768572089794029,
    'negative_predictive_value': 1.0}