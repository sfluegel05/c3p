"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide has repeated units of monosaccharides linked by glycosidic bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for a generic glycosidic linkage (e.g., -O-C-)
    glycosidic_pattern = Chem.MolFromSmarts("O-C")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage found"
    
    # Look for multiple sugar ring units (more flexible pattern)
    cyclic_sugar_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1")  # Simplified monosaccharide pattern
    repeat_units = mol.GetSubstructMatches(cyclic_sugar_pattern)
    if len(repeat_units) < 2:
        return False, "Too few repeating monosaccharide units for polysaccharide"
    
    # Ensure the structure is large enough to be polysaccharide
    if len(mol.GetAtoms()) < 50:  # Arbitrarily chosen, can be adjusted
        return False, "Molecule too small for polysaccharide"
    
    # Polysaccharides should have extensive oxygen atoms due to many OH and glycosidic bonds
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygens < 10:  # Arbitrarily chosen, can be adjusted based on known structures
        return False, "Too few oxygen atoms suggest an incomplete or incorrect polysaccharide"

    return True, "Contains multiple repeated monosaccharide units linked by glycosidic bonds"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18154',
                          'name': 'polysaccharide',
                          'definition': 'A biomacromolecule consisting of '
                                        'large numbers of monosaccharide '
                                        'residues linked glycosidically. This '
                                        'term is commonly used only for those '
                                        'containing more than ten '
                                        'monosaccharide residues.',
                          'parents': [   'CHEBI:16646',
                                         'CHEBI:167559',
                                         'CHEBI:33694'],
                          'xrefs': ['KEGG:C00420'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'O([C@H]1O[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O)[C@@H]1O[C@@H]4O[C@@H]([C@@H](O[C@@H]5O[C@@H]([C@H](O)[C@H](O[C@H]6O[C@@H]([C@H](O)[C@H](O)[C@H]6O)CO)[C@H]5O)CO)[C@H](O)[C@H]4NC(=O)C)CO)CO)[C@H]7[C@H](O)[C@H](O[C@@H](O[C@H]8[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]8CO)O[C@@H]([C@H](O)[C@@H](NC(=O)C)CO)[C@H](O)CO[C@@H]9O[C@H]([C@@H](O)[C@@H](O)[C@@H]9O)C)[C@H]7O)CO[C@H]%10O[C@@H]([C@@H](O)[C@H](O)[C@@H]%10O[C@@H]%11O[C@@H]([C@@H](O[C@@H]%12O[C@@H]([C@H](O)[C@H](O[C@H]%13O[C@@H]([C@H](O)[C@H](O)[C@H]%13O)CO)[C@H]%12O)CO)[C@H](O)[C@H]%11NC(=O)C)CO)CO '
               'NAME: CID 91848939 REASON: MISSED Too few repeating units for '
               'polysaccharide\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)[C@H](O[C@@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@H]3[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]3CO)O)[C@H]1O)CO[C@H]4O[C@@H]([C@@H](O)[C@H](O)[C@@H]4OC5O[C@@H]([C@@H](OC6O[C@@H]([C@H](O)[C@H](OC7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO)[C@H](O)[C@H]5NC(=O)C)CO)CO)[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8OC9O[C@@H]([C@@H](OC%10O[C@@H]([C@H](O)[C@H](OC%11(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%11)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]%10O)CO)[C@H](O)[C@H]9NC(=O)C)CO)CO '
               'NAME: '
               '(4S,5R,6R)-5-Acetamido-2-[(3R,4S,5S,6R)-2-[(2R,3S,4R,5R)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3R,4S,5S,6S)-4-[(2R,3S,4S,5S,6R)-3-[(3R,4R,5S,6R)-3-acetamido-5-[(3R,4S,5S,6R)-4-[(4S,5R,6R)-5-acetamido-2-carboxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxan-2-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,5-dihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
               'acid REASON: MISSED Too few repeating units for '
               'polysaccharide\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@@H]3[C@@H](O[C@H]4[C@H](O)[C@H](O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]5CO)O[C@H]6[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]6CO[C@@H]7O[C@H]([C@@H](O)[C@@H](O)[C@@H]7O)C)O)[C@H]4O)CO[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O[C@@H]%11O[C@@H]([C@@H](O[C@@H]%12O[C@@H]([C@H](O)[C@H](O)[C@H]%12O)CO)[C@H](O[C@@H]%13O[C@H]([C@@H](O)[C@@H](O)[C@@H]%13O)C)[C@H]%11NC(=O)C)CO)[C@H]%10O)CO)[C@H](O[C@@H]%14O[C@H]([C@@H](O)[C@@H](O)[C@@H]%14O)C)[C@H]9NC(=O)C)CO)CO)O[C@@H]([C@@H](O)[C@@H]3O)CO)[C@@H]%15O[C@@H]([C@H](O)[C@H](O[C@@H]%16O[C@@H]([C@@H](O[C@@H]%17O[C@@H]([C@H](O)[C@H](O)[C@H]%17O)CO)[C@H](O[C@@H]%18O[C@H]([C@@H](O)[C@@H](O)[C@@H]%18O)C)[C@H]%16NC(=O)C)CO)[C@H]%15O)CO '
               'NAME: CID 91862369 REASON: MISSED Too few repeating units for '
               'polysaccharide\n'
               ' * SMILES: '
               'O([C@H]1O[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O)[C@@H]1O[C@@H]4O[C@@H]([C@@H](O[C@@H]5O[C@@H]([C@H](O)[C@H](O)[C@H]5O)CO[C@]6(O[C@H]([C@H](NC(=O)C)[C@@H](O)C6)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]4NC(=O)C)CO)CO)[C@H]7[C@H](O)[C@H](O[C@@H](O[C@H]8[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]8CO)O[C@H]9[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]9CO[C@@H]%10O[C@H]([C@@H](O)[C@@H](O)[C@@H]%10O)C)O)[C@H]7O)CO[C@H]%11O[C@@H]([C@@H](O)[C@H](O)[C@@H]%11O[C@@H]%12O[C@@H]([C@@H](O[C@@H]%13O[C@@H]([C@H](O)[C@H](O)[C@H]%13O)CO)[C@H](O)[C@H]%12NC(=O)C)CO)CO[C@@H]%14O[C@@H]([C@@H](O[C@@H]%15O[C@@H]([C@H](O)[C@H](O)[C@H]%15O)CO)[C@H](O)[C@H]%14NC(=O)C)CO '
               'NAME: CID 91857080 REASON: MISSED Too few repeating units for '
               'polysaccharide\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@@H]1CO)OC[C@H]5O[C@@H](O[C@H]6[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]6CO)O[C@@H]([C@H](O)[C@@H](NC(=O)C)CO)[C@H](O)CO)[C@@H](O)[C@@H](O[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@@H]7O[C@@H]8O[C@@H]([C@@H](O[C@@H]9O[C@@H]([C@H](O)[C@H](O)[C@H]9O)CO[C@]%10(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%10)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]8NC(=O)C)CO)CO[C@@H]%11O[C@@H]([C@@H](O[C@@H]%12O[C@@H]([C@H](O)[C@H](O)[C@H]%12O)CO[C@]%13(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%13)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]%11NC(=O)C)CO)[C@@H]5O)[C@@H]%14O[C@@H]([C@@H](O[C@@H]%15O[C@@H]([C@H](O)[C@H](O[C@]%16(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%16)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]%15O)CO)[C@H](O)[C@H]%14NC(=O)C)CO '
               'NAME: CID 91852081 REASON: MISSED Too few repeating units for '
               'polysaccharide\n'
               ' * SMILES: '
               'O1[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@H]3O[C@@H]([C@@H](O)[C@H](O)[C@@H]3O[C@@H]4O[C@@H]([C@@H](O[C@@H]5O[C@@H]([C@H](O)[C@H](O)[C@H]5O)CO)[C@H](O)[C@H]4NC(=O)C)CO)CO)[C@H](O)[C@@H]1O[C@H]6[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]6CO)O[C@@H]([C@H](O)[C@@H](NC(=O)C)CO)[C@H](O)CO[C@@H]7O[C@H]([C@@H](O)[C@@H](O)[C@@H]7O)C)CO[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO)[C@H](O[C@@H]%11O[C@H]([C@@H](O)[C@@H](O)[C@@H]%11O)C)[C@H]9NC(=O)C)CO)CO '
               'NAME: '
               'N-[(2S,3R,4R,5S,6R)-2-[(2R,3R,4R,5S,6S)-6-[(2R,3S,4R,5R,6S)-5-Acetamido-6-[(2R,3S,4R,5S)-5-acetamido-2,4,6-trihydroxy-1-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyhexan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-4-[(2R,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-4-hydroxy-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2-[[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-5-hydroxyoxan-3-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide '
               'REASON: MISSED Too few repeating units for polysaccharide\n'
               ' * SMILES: '
               '[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O[C@@H]2[C@H]([C@@H](O[C@@H]([C@H]2O)CO)O[C@@H]3[C@H]([C@@H](O[C@@H]([C@H]3O)CO)O[C@@H]4[C@H]([C@@H](O[C@@H]([C@H]4O)CO)O[C@@H]5[C@H]([C@@H](O[C@@H]([C@H]5O)CO)O[C@@H]6[C@H]([C@@H](O[C@@H]([C@H]6O)CO)O[C@@H]7[C@H]([C@@H](O[C@@H]([C@H]7O)CO)O[C@@H]8[C@H]([C@@H](O[C@@H]([C@H]8O)CO)O[C@@H]9[C@H]([C@@H](O[C@@H]([C@H]9O)CO)O[C@@H]%10[C@H]([C@@H](O[C@@H]([C@H]%10O)CO)O[C@@H]%11[C@H]([C@@H](O[C@@H]([C@H]%11O)CO)O[C@@H]%12[C@H]([C@@H](O[C@@H]([C@H]%12O)CO)O)O)O)O)O)O)O)O)O)O)O)O '
               'NAME: (1->3)-dodeca-beta-D-glucan REASON: MISSED Too few '
               'repeating units for polysaccharide\n'
               ' * SMILES: '
               'O([C@@H]1[C@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@H]3[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]3CO[C@H]4O[C@@H]([C@H](O)[C@H](O)[C@H]4O)C)O)O[C@@H]([C@H]1O)CO[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O)CO[C@]8(O[C@H]([C@H](NC(=O)C)[C@@H](O)C8)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]6NC(=O)C)CO)CO)[C@H]9O[C@@H]([C@@H](O)[C@H](O)[C@@H]9O[C@@H]%10O[C@@H]([C@@H](O[C@@H]%11O[C@@H]([C@H](O)[C@H](O)[C@H]%11O)CO[C@]%12(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%12)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]%10NC(=O)C)CO)CO '
               'NAME: '
               '(2R,4S,5R,6R)-5-Acetamido-2-[[(2R,3R,4S,5R,6S)-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3R,4S,5S,6S)-4-[(2R,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-5-[(2S,3R,4S,5R,6R)-6-[[(2R,4S,5R,6R)-5-acetamido-2-carboxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxan-2-yl]oxymethyl]-3,4,5-trihydroxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5R,6R)-5-acetamido-4,6-dihydroxy-2-[[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxymethyl]oxan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,5-dihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
               'acid REASON: MISSED Too few repeating units for '
               'polysaccharide\n'
               ' * SMILES: '
               'O1[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@H]3O[C@@H]([C@@H](O)[C@H](O)[C@@H]3O[C@@H]4O[C@@H]([C@@H](O[C@@H]5O[C@@H]([C@H](O)[C@H](O)[C@H]5O)CO[C@]6(O[C@H]([C@H](NC(=O)C)[C@@H](O)C6)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]4NC(=O)C)CO)CO)[C@H](O)[C@@H]1O[C@H]7[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]7CO)O[C@H]8[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]8CO[C@@H]9O[C@H]([C@@H](O)[C@@H](O)[C@@H]9O)C)O)CO[C@H]%10O[C@@H]([C@@H](O)[C@H](O)[C@@H]%10O[C@@H]%11O[C@@H]([C@@H](O[C@@H]%12O[C@@H]([C@H](O)[C@H](O[C@]%13(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%13)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]%12O)CO)[C@H](O)[C@H]%11NC(=O)C)CO)CO '
               'NAME: CID 71297742 REASON: MISSED Too few repeating units for '
               'polysaccharide\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)O)[C@@H]4O[C@@H]([C@@H](O)[C@H](O[C@H]5O[C@@H]([C@@H](O[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O)CO)[C@H](O)[C@H]6NC(=O)C)CO)[C@H](O)[C@@H]5O)CO)[C@@H]4O)CO[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO)[C@H](O)[C@H]9NC(=O)C)CO)CO[C@@H]%11O[C@@H]([C@@H](O[C@@H]%12O[C@@H]([C@H](O)[C@H](O)[C@H]%12O)CO)[C@H](O)[C@H]%11NC(=O)C)CO '
               'NAME: '
               'N-[(2R,3R,4R,5S,6R)-5-[(2S,3R,4R,5S,6R)-3-Acetamido-5-[(2S,3S,4S,5R,6R)-6-[[(2S,3S,4S,5S,6R)-3-[(2R,3R,4R,5S,6R)-3-acetamido-4-hydroxy-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-6-[[(2R,3R,4R,5S,6R)-3-acetamido-4-hydroxy-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-4,5-dihydroxyoxan-2-yl]oxymethyl]-4-[(2R,3S,4R,5S,6R)-5-[(2R,3R,4R,5S,6R)-3-acetamido-4-hydroxy-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-3,4-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-3,5-dihydroxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-[[(2R,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxymethyl]oxan-3-yl]acetamide '
               'REASON: MISSED Too few repeating units for polysaccharide\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'Too few repeating '
                                               'monosaccharide units for '
                                               'polysaccharide'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'Too few repeating '
                                               'monosaccharide units for '
                                               'polysaccharide'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'No glycosidic linkage found'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'Too few repeating '
                                               'monosaccharide units for '
                                               'polysaccharide'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'No glycosidic linkage found'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'Too few repeating '
                                               'monosaccharide units for '
                                               'polysaccharide'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'Too few repeating '
                                               'monosaccharide units for '
                                               'polysaccharide'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'Too few repeating '
                                               'monosaccharide units for '
                                               'polysaccharide'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'Too few repeating '
                                               'monosaccharide units for '
                                               'polysaccharide'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'Molecule too small for '
                                               'polysaccharide'}],
    'sample_false_negatives': [   {   'smiles': 'O(C1=C(C(OC([H])([H])[H])=C(C(=C1[H])/C(=C(/C2=C(C(O[H])=C(O[H])C(=C2[H])[H])[H])\\[H])/[H])[H])\\C(=C(\\C(C([H])([H])[H])(C([H])([H])[H])[H])/[H])\\[H])C([H])([H])[H]',
                                      'name': 'mucilages A',
                                      'reason': 'Too few repeating '
                                                'monosaccharide units for '
                                                'polysaccharide'}],
    'attempt': 1,
    'reasoning': 'The previous attempt to classify a polysaccharide was likely '
                 'ineffective due to several issues that need to be refined. '
                 "Here's what could have possibly gone wrong and how to "
                 'address it:\n'
                 '\n'
                 '1. **Glycosidic Linkage Identification**: The glycosidic '
                 'bond pattern that was used may not be specific or '
                 'comprehensive enough to capture a variety of glycosidic '
                 'linkages. Glycosidic bonds can vary significantly depending '
                 'on the type of sugar and the linkage (α or β).\n'
                 '\n'
                 '2. **Cyclic Monosaccharide Pattern**: The pattern used to '
                 'identify cyclic monosaccharide units might have been too '
                 'specific. Polysaccharides can consist of a variety of '
                 'monosaccharides, and slight variations in their '
                 'stereochemistry might have prevented accurate matching.\n'
                 '\n'
                 '3. **Number of Units as Criterion**: The criterion "too few '
                 'repeating units" might need additional adjustments to handle '
                 'larger or differently structured polysaccharides. We could '
                 'consider more flexible criteria for counting units depending '
                 'on the compound or relax the constraints if needed.\n'
                 '\n'
                 '4. **Hydroxyl Groups Detection**: Simply counting the '
                 'hydroxyl groups may not account for the complexity and '
                 'variability in polysaccharide structures. \n'
                 '\n'
                 '5. **Molecular Weight Check**: Polysaccharides can vary '
                 'greatly in size; thus, using a strict cutoff based on '
                 'molecular weight alone may exclude valid polysaccharide '
                 'structures.\n'
                 '\n'
                 'To improve the classification, we can:\n'
                 '\n'
                 '- Broaden the definition of glycosidic bonds using more '
                 'comprehensive SMARTS patterns that capture the typical -O-C- '
                 'linkage between sugars, accounting for various anomers.\n'
                 '- Capture the cyclic backbone structure of polysaccharides '
                 'through more general patterns that recognize carbohydrate '
                 'rings without overly strict stereochemistry constraints.\n'
                 '- Redefine the "repeating units" requirement to account for '
                 'the diversity seen across different polysaccharides.\n'
                 '- Rather than relying solely on hydroxyl counts or molecular '
                 'weight, we can focus on connectivity patterns between sugar '
                 'units.\n'
                 '\n'
                 "Here's an improved version of the code considering these "
                 'factors:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 1296,
    'num_false_positives': 100,
    'num_true_negatives': 2152,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.9283667621776505,
    'recall': 0.9992289899768697,
    'f1': 0.9624953583364277,
    'accuracy': 0.9715412792335869,
    'negative_predictive_value': 0.9995355318160706}