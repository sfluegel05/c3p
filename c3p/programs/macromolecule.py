"""
Classifies: CHEBI:33839 macromolecule
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    A macromolecule is defined as having a high molecular mass,
    and the structure comprised of repetitive units from smaller molecules.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a macromolecule, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Check for high molecular weight - typically in 5000 Da or greater
    if mol_wt >= 5000:
        return True, f"Molecule has high molecular weight: {mol_wt:.2f} Da"
    
    # Count total number of atoms
    atom_count = mol.GetNumHeavyAtoms()
    
    # Check for repetition patterns, using common biological motifs like peptide linkage
    # Implement checks for repeating units such as -[CH2-]n, amide bonds for peptides, etc.
    peptide_pattern = Chem.MolFromSmarts("N[C](=O)[C]")
    is_peptide_repetitive = len(mol.GetSubstructMatches(peptide_pattern)) > 3
    
    # Approximate check for a large number of atoms and repeating subunits might denote a polymeric nature
    if atom_count >= 100 or is_peptide_repetitive:
        return True, f"Molecule has a large number of atoms ({atom_count}) or peptide-like repeating units"
    
    return False, "Molecule does not meet criteria for a macromolecule"

# Example testing
example_smiles = "O1[C@H](O[C@H]2OC(CO)C(O)C2O[C@@H](C(O)CO)O[C@@H]3O[C@@H]([C@@H](O)C(O)[C@H]3O)CO)C(O)C(O)[C@H](CO)O[C@H]4OC(CO)C(O)C4O"
print(is_macromolecule(example_smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33839',
                          'name': 'macromolecule',
                          'definition': 'A macromolecule is a molecule of high '
                                        'relative molecular mass, the '
                                        'structure of which essentially '
                                        'comprises the multiple repetition of '
                                        'units derived, actually or '
                                        'conceptually, from molecules of low '
                                        'relative molecular mass.',
                          'parents': ['CHEBI:36357'],
                          'xrefs': ['Wikipedia:Macromolecule'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.8975409836065573 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C '
               'NAME: Chrysospermin B REASON: CORRECT Molecule has a large '
               'number of atoms (136) or repeating substructures (4)\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)[C@H](O[C@@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@H]3[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]3CO[C@@H]4O[C@H]([C@@H](O)[C@@H](O)[C@@H]4O)C)O)[C@H]1O)CO[C@H]5O[C@@H]([C@@H](O[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O[C@@H]8O[C@@H]([C@H](O)[C@H](O)[C@H]8O)CO)[C@H]7O)CO)[C@H](O)[C@H]6NC(=O)C)CO)[C@H](O)[C@@H]5O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO)[C@H](O)[C@H]9NC(=O)C)CO)CO)[C@H]%11O[C@@H]([C@@H](O)[C@H](O)[C@@H]%11O[C@@H]%12O[C@@H]([C@@H](O[C@@H]%13O[C@@H]([C@H](O)[C@H](O)[C@H]%13O)CO)[C@H](O)[C@H]%12NC(=O)C)CO)CO[C@@H]%14O[C@@H]([C@@H](O[C@@H]%15O[C@@H]([C@H](O)[C@H](O)[C@H]%15O)CO)[C@H](O)[C@H]%14NC(=O)C)CO '
               'NAME: CID 91846516 REASON: CORRECT Molecule has a large number '
               'of atoms (183) or repeating substructures (15)\n'
               ' * SMILES: '
               'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC2=CC=CC=C2)CCC(=O)N)CCC(=O)N)(CC)C)(C)C)C(C)C)CCC1)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)(C)C)C)(C)C)C)(C)C)(C)C)CCC(=O)N)(C)C)C(C)C)(C)C)CC(C)C)(C)C '
               'NAME: Trichokonin-IX REASON: CORRECT Molecule has a large '
               'number of atoms (139) or repeating substructures (2)\n'
               ' * SMILES: '
               'O=C1NC(C(=O)NCC(=O)NC(C(=O)NC2C(=O)NC(C(=O)NC(CC=3N=CNC3)C(NC(C(NC(C(NC(CC(OCC(C(NC1C(O)C)=O)NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(N)CCC(=O)O)CC(C)C)C(O)C)C(O)C)C(O)C)=O)C(=O)O)=O)CC(=O)NCCCC2)=O)CC4=CC=CC=C4)=O)C(O)C)CO)C '
               'NAME: Prunipeptin REASON: CORRECT Molecule has a large number '
               'of atoms (122) or repeating substructures (4)\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)[C@H](O[C@@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@H]3[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]3CO[C@@H]4O[C@H]([C@@H](O)[C@@H](O)[C@@H]4O)C)O)[C@H]1O)CO[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O[C@@H]8O[C@@H]([C@@H](O[C@@H]9O[C@@H]([C@H](O)[C@H](O)[C@H]9O)CO)[C@H](O)[C@H]8NC(=O)C)CO)[C@H]7O)CO)[C@H](O)[C@H]6NC(=O)C)CO)CO[C@@H]%10O[C@@H]([C@@H](O[C@@H]%11O[C@@H]([C@H](O)[C@H](O)[C@H]%11O)CO)[C@H](O)[C@H]%10NC(=O)C)CO)[C@H]%12O[C@@H]([C@@H](O[C@@H]%13O[C@@H]([C@@H](O[C@@H]%14O[C@@H]([C@H](O)[C@H](O)[C@H]%14O)CO)[C@H](O)[C@H]%13NC(=O)C)CO)[C@H](O)[C@@H]%12O[C@@H]%15O[C@@H]([C@@H](O[C@@H]%16O[C@@H]([C@H](O)[C@H](O)[C@H]%16O)CO)[C@H](O)[C@H]%15NC(=O)C)CO)CO '
               'NAME: CID 91853261 REASON: CORRECT Molecule has a large number '
               'of atoms (197) or repeating substructures (16)\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@@H]1CO)O[C@H]5[C@H](O)[C@H](O[C@@H](O[C@H]6[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]6CO)O[C@@H]([C@H](O)[C@@H](NC(=O)C)CO)[C@H](O)CO)[C@H]5O)CO[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@@H]7O[C@@H]8O[C@@H]([C@@H](O[C@@H]9O[C@@H]([C@H](O)[C@H](O)[C@H]9O)CO[C@]%10(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%10)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]8NC(=O)C)CO)CO)[C@@H]%11O[C@@H]([C@@H](O)[C@H](O[C@@H]%12O[C@@H]([C@H](O)[C@H](O)[C@H]%12O)CO)[C@H]%11NC(=O)C)CO[C@]%13(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%13)[C@H](O)[C@H](O)CO)C(O)=O '
               'NAME: CID 91859078 REASON: CORRECT Molecule has a large number '
               'of atoms (197) or repeating substructures (13)\n'
               ' * SMILES: '
               'O([C@@H]1[C@H](O)[C@@H](O[C@@H]([C@H]1O)CO[C@H]2O[C@@H]([C@@H](O)[C@H](O)[C@@H]2O[C@@H]3O[C@@H]([C@@H](O[C@@H]4O[C@@H]([C@H](O)[C@H](O)[C@H]4O[C@@H]5O[C@H]([C@@H](O)[C@@H](O)[C@@H]5O)C)CO)[C@H](O)[C@H]3NC(=O)C)CO)CO)O[C@H]6[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]6CO)O[C@H]7[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]7CO)O)[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O[C@@H]%11O[C@H]([C@@H](O)[C@@H](O)[C@@H]%11O)C)CO)[C@H](O)[C@H]9NC(=O)C)CO)CO '
               'NAME: '
               'N-[(3R,4R,5S,6R)-5-[(2S,3R,4R,5S,6R)-3-Acetamido-5-[(2S,3S,4S,5R,6R)-4-[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-5-[(2S,3R,4S,5R,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-5-[(2S,3R,4S,5R,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-3,5-dihydroxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide '
               'REASON: CORRECT Molecule has a large number of atoms (132) or '
               'repeating substructures (11)\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@@H](O[C@H]3[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]3CO)O[C@H]4[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]4CO)O)[C@H]1O)CO[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O[C@@H]8O[C@@H]([C@@H](O[C@@H]9O[C@@H]([C@H](O)[C@H](O[C@@H]%10O[C@@H]([C@@H](O[C@@H]%11O[C@@H]([C@H](O)[C@H](O[C@@H]%12O[C@@H]([C@@H](O[C@@H]%13O[C@@H]([C@H](O)[C@H](O[C@@H]%14O[C@@H]([C@@H](O[C@@H]%15O[C@@H]([C@H](O)[C@H](O[C@@H]%16O[C@@H]([C@@H](O[C@@H]%17O[C@@H]([C@H](O)[C@H](O)[C@H]%17O)CO)[C@H](O)[C@H]%16NC(=O)C)CO)[C@H]%15O)CO)[C@H](O)[C@H]%14NC(=O)C)CO)[C@H]%13O)CO)[C@H](O)[C@H]%12NC(=O)C)CO)[C@H]%11O)CO)[C@H](O)[C@H]%10NC(=O)C)CO)[C@H]9O)CO)[C@H](O)[C@H]8NC(=O)C)CO)[C@H]7O)CO)[C@H](O)[C@H]6NC(=O)C)CO)CO)[C@H]%18O[C@@H]([C@@H](O)[C@H](O)[C@@H]%18O[C@@H]%19O[C@@H]([C@@H](O[C@@H]%20O[C@@H]([C@H](O)[C@H](O[C@@H]%21O[C@@H]([C@@H](O[C@@H]%22O[C@@H]([C@H](O)[C@H](O)[C@H]%22O)CO[C@]%23(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%23)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]%21NC(=O)C)CO)[C@H]%20O)CO)[C@H](O)[C@H]%19NC(=O)C)CO)CO '
               'NAME: CID 91847754 REASON: CORRECT Molecule has a large number '
               'of atoms (296) or repeating substructures (23)\n'
               ' * SMILES: '
               'O1[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@H]3O[C@@H]([C@@H](O)[C@H](O)[C@@H]3O[C@@H]4O[C@@H]([C@@H](O[C@@H]5O[C@@H]([C@H](O)[C@H](O)[C@H]5O)CO)[C@H](O)[C@H]4NC(=O)C)CO)CO)[C@H](O)[C@@H]1O[C@H]6[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]6CO)O[C@H]7[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]7CO)O)CO[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO)[C@H](O[C@@H]%11O[C@H]([C@@H](O)[C@@H](O)[C@@H]%11O)C)[C@H]9NC(=O)C)CO)CO '
               'NAME: '
               'N-[(2R,3R,4R,5S,6R)-5-[(2S,3R,4R,5S,6R)-3-Acetamido-5-[(2S,3S,4R,5R,6R)-5-[(2S,3R,4R,5S,6R)-3-acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-[(2R,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-4-hydroxy-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-3-hydroxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide '
               'REASON: CORRECT Molecule has a large number of atoms (136) or '
               'repeating substructures (11)\n'
               ' * SMILES: '
               'O([C@H]1[C@@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@@H]1OC[C@H]5O[C@@H](O[C@H]6[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]6CO)O[C@H]7[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]7CO)O)[C@@H](O)[C@@H](O[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO)[C@H](O)[C@H]9NC(=O)C)CO)CO)[C@@H]5O)CO)[C@@H]%11O[C@@H]([C@@H](O[C@@H]%12O[C@@H]([C@H](O)[C@H](O)[C@H]%12O)CO[C@]%13(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%13)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]%11NC(=O)C)CO '
               'NAME: CID 91857085 REASON: CORRECT Molecule has a large number '
               'of atoms (177) or repeating substructures (13)\n'
               ' * SMILES: '
               'O([C@H]1O[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]3O)CO)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O)[C@@H]1O[C@@H]5O[C@@H]([C@@H](O[C@@H]6O[C@@H]([C@H](O)[C@H](O)[C@H]6O)CO)[C@H](O)[C@H]5NC(=O)C)CO)CO)[C@H]7[C@H](O)[C@H](O[C@@H](O[C@H]8[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]8CO)O[C@H]9[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]9CO[C@@H]%10O[C@H]([C@@H](O)[C@@H](O)[C@@H]%10O)C)O)[C@H]7O)CO[C@H]%11O[C@@H]([C@@H](O)[C@H](O)[C@@H]%11O[C@@H]%12O[C@@H]([C@@H](O[C@@H]%13O[C@@H]([C@H](O)[C@H](O[C@]%14(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%14)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]%13O)CO)[C@H](O)[C@H]%12NC(=O)C)CO)CO[C@@H]%15O[C@@H]([C@@H](O[C@@H]%16O[C@@H]([C@H](O)[C@H](O[C@]%17(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%17)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]%16O)CO)[C@H](O)[C@H]%15NC(=O)C)CO '
               'NAME: CID 91850256 REASON: CORRECT Molecule has a large number '
               'of atoms (232) or repeating substructures (17)\n'
               ' * SMILES: '
               'O([C@@H]1[C@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@H]3[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]3CO[C@@H]4O[C@H]([C@@H](O)[C@@H](O)[C@@H]4O)C)O)O[C@@H]([C@H]1O)CO[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O)CO)[C@H](O)[C@H]6NC(=O)C)CO)CO[C@@H]8O[C@@H]([C@@H](O[C@@H]9O[C@@H]([C@H](O)[C@H](O)[C@H]9NC(=O)C)CO)[C@H](O)[C@H]8NC(=O)C)CO)[C@H]%10O[C@@H]([C@@H](O)[C@H](O)[C@@H]%10O[C@@H]%11O[C@@H]([C@@H](O[C@@H]%12O[C@@H]([C@H](O)[C@H](O)[C@H]%12NC(=O)C)CO)[C@H](O)[C@H]%11NC(=O)C)CO)CO '
               'NAME: '
               'N-[(2R,3R,4R,5S,6R)-5-[(2S,3R,4R,5S,6R)-3-Acetamido-5-[(2S,3S,4S,5R,6R)-4-[(2R,3S,4S,5S,6R)-3-[(2R,3R,4R,5S,6R)-3-acetamido-5-[(2S,3R,4R,5R,6R)-3-acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2S,3S,4S,5S,6R)-6-[[(2R,3R,4R,5S,6R)-3-acetamido-5-[(2S,3R,4R,5R,6R)-3-acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-3-[(2S,3R,4R,5S,6R)-3-acetamido-4-hydroxy-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-4,5-dihydroxyoxan-2-yl]oxymethyl]-3,5-dihydroxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-[[(2R,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxymethyl]oxan-3-yl]acetamide '
               'REASON: CORRECT Molecule has a large number of atoms (153) or '
               'repeating substructures (12)\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]2[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]2CO[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)O)[C@@H]4O[C@@H]([C@@H](O)[C@H](O[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O)CO)[C@H](O)[C@H]6NC(=O)C)CO)CO)[C@@H]4O)CO[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O[C@@H]%11O[C@@H]([C@@H](O[C@@H]%12O[C@@H]([C@H](O)[C@H](O[C@@H]%13O[C@@H]([C@@H](O[C@@H]%14O[C@@H]([C@H](O)[C@H](O)[C@H]%14O)CO)[C@H](O)[C@H]%13NC(=O)C)CO)[C@H]%12O)CO)[C@H](O)[C@H]%11NC(=O)C)CO)[C@H]%10O)CO)[C@H](O)[C@H]9NC(=O)C)CO)CO '
               'NAME: CID 91855423 REASON: CORRECT Molecule has a large number '
               'of atoms (172) or repeating substructures (14)\n'
               ' * SMILES: '
               'O([C@H]1[C@@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]3O)CO)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@@H]1O[C@H]5[C@H](O)[C@H](O[C@@H](O[C@H]6[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]6CO)O[C@H]7[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]7CO)O)[C@H]5O)CO[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO[C@]%11(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%11)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]9NC(=O)C)CO)CO)CO)[C@@H]%12O[C@@H]([C@@H](O[C@@H]%13O[C@@H]([C@H](O)[C@H](O[C@]%14(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%14)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]%13O)CO)[C@H](O)[C@H]%12NC(=O)C)CO '
               'NAME: '
               '(2R,4S,5R,6R)-5-Acetamido-2-[[(2R,3R,4S,5R,6S)-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3R,4S,5S,6S)-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-4-[(2R,3S,4S,5S,6R)-3,5-bis[[(2S,3R,4R,5S,6R)-3-acetamido-5-[(2S,3R,4S,5S,6R)-4-[(2S,4S,5R,6R)-5-acetamido-2-carboxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxan-2-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy]-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-3,5-dihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
               'acid REASON: CORRECT Molecule has a large number of atoms '
               '(197) or repeating substructures (14)\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]2[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]2CO)O)[C@H]3O[C@@H]([C@@H](O)[C@H](O[C@H]4O[C@@H]([C@@H](O)[C@H](O)[C@@H]4O[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O[C@H]6O[C@@H]([C@@H](O)[C@H](O)[C@@H]6O)CO)CO)CO)[C@@H]3O)CO[C@H]7O[C@@H]([C@@H](O)[C@H](O[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@H]9O[C@@H]([C@@H](O)[C@H](O)[C@@H]9O)CO)CO)[C@@H]7O)CO[C@H]%10O[C@@H]([C@@H](O)[C@H](O)[C@@H]%10O[C@H]%11O[C@@H]([C@@H](O)[C@H](O)[C@@H]%11O)CO)CO '
               'NAME: '
               'N-[(3R,4R,5S,6R)-5-[(2S,3R,4R,5S,6R)-3-Acetamido-5-[(2R,3S,4S,5R,6R)-4-[(2R,3S,4S,5S,6R)-3-[(2S,3S,4S,5S,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2R,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2S,3S,4S,5R,6R)-4-[(2S,3S,4S,5S,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2R,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-6-[[(2S,3S,4S,5S,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2R,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-3,5-dihydroxyoxan-2-yl]oxymethyl]-3,5-dihydroxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide '
               'REASON: CORRECT Molecule has a large number of atoms (128) or '
               'repeating substructures (11)\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@@H](O[C@H]3[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]3CO)O[C@H]4[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]4CO)O)[C@H]1O)CO[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O)CO[C@]8(O[C@H]([C@H](NC(=O)C)[C@@H](O)C8)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]6NC(=O)C)CO)CO)[C@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@@H](O[C@@H]%11O[C@@H]([C@H](O[C@H]%12O[C@@H]([C@H](O)[C@H](O)[C@H]%12O)CO)[C@H](O)[C@H]%11O)CO)[C@H](O)[C@H]%10NC(=O)C)CO)[C@H](O)[C@@H]9O[C@@H]%13O[C@@H]([C@@H](O[C@@H]%14O[C@@H]([C@H](O)[C@H](O)[C@H]%14O)CO[C@]%15(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%15)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]%13NC(=O)C)CO)CO '
               'NAME: CID 71297944 REASON: CORRECT Molecule has a large number '
               'of atoms (202) or repeating substructures (15)\n'
               ' * SMILES: '
               'O([C@@H]1[C@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@H]3[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]3CO)O)O[C@@H]([C@H]1O)CO[C@H]4O[C@@H]([C@@H](O)[C@H](O[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O)CO)[C@@H]4O)CO[C@H]6O[C@@H]([C@@H](O)[C@H](O)[C@@H]6O)CO)[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@@H]7O[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@H]9O[C@@H]([C@@H](O)[C@H](O[C@H]%10O[C@@H]([C@@H](O)[C@H](O)[C@@H]%10O)CO)[C@@H]9O)CO)CO)CO[C@H]%11O[C@@H]([C@@H](O)[C@H](O)[C@@H]%11O[C@H]%12O[C@@H]([C@@H](O)[C@H](O)[C@@H]%12O)CO)CO '
               'NAME: '
               'N-[(3R,4R,5S,6R)-5-[(2S,3R,4R,5S,6R)-3-Acetamido-5-[(2S,3S,4S,5R,6R)-4-[(2R,3S,4S,5S,6R)-3-[(2R,3S,4S,5S,6R)-3-[(2R,3S,4S,5R,6R)-3,5-dihydroxy-6-(hydroxymethyl)-4-[(2S,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2S,3S,4S,5S,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2S,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-4,5-dihydroxyoxan-2-yl]oxy-6-[[(2S,3S,4S,5R,6R)-3,5-dihydroxy-4-[(2S,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2S,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-2-yl]oxymethyl]-3,5-dihydroxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide '
               'REASON: CORRECT Molecule has a large number of atoms (139) or '
               'repeating substructures (12)\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@@H]1CO)O[C@H]4[C@H](O)[C@H](O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]5CO)O[C@H]6[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]6CO)O)[C@H]4O)CO[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@@H]7O[C@@H]8O[C@@H]([C@@H](O[C@@H]9O[C@@H]([C@H](O)[C@H](O[C@]%10(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%10)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]9O)CO)[C@H](O)[C@H]8NC(=O)C)CO)CO)[C@@H]%11O[C@@H]([C@@H](O)[C@H](O[C@@H]%12O[C@@H]([C@H](O)[C@H](O[C@]%13(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%13)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]%12O)CO)[C@H]%11NC(=O)C)CO[C@]%14(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%14)[C@H](O)[C@H](O)CO)C(O)=O '
               'NAME: CID 91862275 REASON: CORRECT Molecule has a large number '
               'of atoms (197) or repeating substructures (14)\n'
               ' * SMILES: '
               'O([C@@H]1[C@@H](O)[C@H](O[C@@H]2[C@@H](NC(=O)C)[C@@H](O[C@@H]([C@H]2O)CO)O[C@@H]3[C@@H](O)[C@@H](O[C@@H]([C@@H]3O)CO)O[C@@H]4[C@@H](NC(=O)C)C(O[C@@H]([C@@H]4O)CO)O)O[C@@H]([C@@H]1O)CO[C@@H]5O[C@@H]([C@@H](O)[C@H](O[C@@H]6O[C@@H]([C@H](O)[C@H](O)[C@H]6O[C@@H]7O[C@H]([C@@H](O)[C@@H](O)[C@@H]7O)C)CO)[C@H]5NC(=O)C)CO)[C@@H]8O[C@@H]([C@@H](O)[C@H](O[C@@H]9O[C@@H]([C@H](O)[C@H](O[C@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10NC(=O)C)CO)[C@H]9O[C@@H]%11O[C@H]([C@@H](O)[C@@H](O)[C@@H]%11O)C)CO)[C@H]8NC(=O)C)CO '
               'NAME: '
               'N-[(2R,3R,4R,5R,6R)-2-[(2R,3R,4S,5S,6R)-2-[(2S,3R,4R,5S,6R)-3-Acetamido-2-[(2R,3R,4S,5S,6R)-2-[(2S,3R,4R,5S,6R)-3-acetamido-2-[(2R,3R,4S,5S,6R)-2-[(3R,4R,5R,6R)-3-acetamido-2,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-5-hydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-6-[[(2R,3R,4R,5S,6R)-3-acetamido-4-[(2R,3R,4S,5R,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-5-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-3,5-dihydroxyoxan-4-yl]oxy-5-hydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-5-hydroxy-6-(hydroxymethyl)-3-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-4-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide '
               'REASON: CORRECT Molecule has a large number of atoms (135) or '
               'repeating substructures (11)\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO)[C@H](O)[C@H]2NC(=O)C)CO)[C@@H](O[C@H]4[C@H](O)[C@H](O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]5CO)O[C@H]6[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]6CO)O)[C@H]4O)CO[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@@H]7O[C@@H]8O[C@@H]([C@@H](O[C@@H]9O[C@@H]([C@H](O)[C@H](O)[C@H]9O)CO)[C@H](O[C@@H]%10O[C@H]([C@@H](O)[C@@H](O)[C@@H]%10O)C)[C@H]8NC(=O)C)CO)CO[C@@H]%11O[C@@H]([C@@H](O[C@@H]%12O[C@@H]([C@H](O)[C@H](O)[C@H]%12O)CO)[C@H](O)[C@H]%11NC(=O)C)CO)O[C@@H]1CO)[C@@H]%13O[C@@H]([C@@H](O[C@@H]%14O[C@@H]([C@H](O)[C@H](O)[C@H]%14O)CO)[C@H](O[C@@H]%15O[C@H]([C@@H](O)[C@@H](O)[C@@H]%15O)C)[C@H]%13NC(=O)C)CO '
               'NAME: CID 91850383 REASON: CORRECT Molecule has a large number '
               'of atoms (182) or repeating substructures (15)\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@@H](O[C@H]3[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]3CO)O[C@H]4[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]4CO[C@@H]5O[C@H]([C@@H](O)[C@@H](O)[C@@H]5O)C)O)[C@H]1O)CO[C@H]6O[C@@H]([C@@H](O)[C@H](O[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@@H]7O)CO)[C@@H]6O)CO[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O)CO)[C@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@@H](O[C@@H]%11O[C@@H]([C@H](O)[C@H](O)[C@H]%11O)CO)[C@H](O)[C@H]%10NC(=O)C)CO)[C@H](O)[C@@H]9O[C@@H]%12O[C@@H]([C@@H](O[C@@H]%13O[C@@H]([C@H](O)[C@H](O)[C@H]%13O)CO)[C@H](O)[C@H]%12NC(=O)C)CO)CO '
               'NAME: '
               'N-[(3R,4R,5S,6R)-5-[(2S,3R,4R,5S,6R)-3-Acetamido-5-[(2S,3S,4R,5R,6R)-5-[(2S,3R,4R,5S,6R)-3-acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-[(2R,3S,4S,5S,6R)-3,5-bis[[(2S,3R,4R,5S,6R)-3-acetamido-4-hydroxy-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy]-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2S,3S,4S,5R,6R)-3,5-dihydroxy-4-[(2S,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2S,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-2-yl]oxymethyl]-3-hydroxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-[[(2R,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxymethyl]oxan-3-yl]acetamide '
               'REASON: CORRECT Molecule has a large number of atoms (158) or '
               'repeating substructures (13)\n'
               ' * SMILES: '
               'O([C@@H]1[C@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@H]3[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]3CO[C@@H]4O[C@H]([C@@H](O)[C@@H](O)[C@@H]4O)C)O)O[C@@H]([C@H]1O)CO[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7NC(=O)C)CO)[C@H](O)[C@H]6NC(=O)C)CO)CO)[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10NC(=O)C)CO[C@]%11(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%11)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]9NC(=O)C)CO)CO '
               'NAME: '
               '(2R,4S,5R,6R)-5-Acetamido-2-[[(2R,3R,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4S,5S,6R)-2-[(2R,3R,4S,5S,6S)-2-[[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-5-[(2S,3R,4R,5R,6R)-3-acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5R,6R)-5-acetamido-4,6-dihydroxy-2-[[(2R,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxymethyl]oxan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,5-dihydroxyoxan-4-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4-dihydroxyoxan-2-yl]methoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
               'acid REASON: CORRECT Molecule has a large number of atoms '
               '(148) or repeating substructures (11)\n'
               ' * SMILES: '
               'O([C@H]1[C@@H](O)[C@H](O[C@@H](O[C@H]2[C@@H](O)[C@H](OC(O)[C@@H]2NC(=O)C)CO[C@@H]3O[C@@H]([C@@H](O[C@@H]4O[C@@H]([C@H](O)[C@H](O)[C@H]4O[C@@H]5O[C@H]([C@@H](O)[C@@H](O)[C@@H]5O)C)CO)[C@H](O)[C@H]3NC(=O)C)CO)[C@@H]1O)CO[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O[C@@H]8O[C@H]([C@@H](O)[C@@H](O)[C@@H]8O)C)CO)[C@H](O)[C@H]6NC(=O)C)CO)[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O[C@@H]%11O[C@H]([C@@H](O)[C@@H](O)[C@@H]%11O)C)CO)[C@H](O)[C@H]9NC(=O)C)CO '
               'NAME: '
               'N-[(2R,3R,4R,5S,6R)-2-[[(2R,3R,4R,5R)-5-Acetamido-4-[(2R,3R,4S,5S,6R)-4-[(2R,3R,4R,5S,6R)-3-acetamido-5-[(2S,3R,4S,5R,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3R,4R,5S,6R)-3-acetamido-5-[(2S,3R,4S,5R,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-3,5-dihydroxyoxan-2-yl]oxy-3,6-dihydroxyoxan-2-yl]methoxy]-5-[(2S,3R,4S,5R,6R)-4,5-dihydroxy-6-(hydroxymethyl)-3-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide '
               'REASON: CORRECT Molecule has a large number of atoms (131) or '
               'repeating substructures (11)\n'
               ' * SMILES: '
               'O([C@@H]1[C@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@H]3[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]3CO[C@@H]4O[C@H]([C@@H](O)[C@@H](O)[C@@H]4O)C)O)O[C@@H]([C@H]1O)CO[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O)CO[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O)CO)[C@H](O)[C@H]6NC(=O)C)CO)[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO)[C@H](O)[C@H]9NC(=O)C)CO)CO[C@@H]%11O[C@@H]([C@@H](O[C@@H]%12O[C@@H]([C@H](O)[C@H](O)[C@H]%12O)CO)[C@H](O)[C@H]%11NC(=O)C)CO '
               'NAME: '
               'N-[(2R,3R,4R,5S,6R)-5-[(2S,3R,4R,5S,6R)-3-Acetamido-5-[(2S,3S,4S,5R,6R)-4-[(2R,3S,4S,5S,6R)-3-[(2R,3R,4R,5S,6R)-3-acetamido-4-hydroxy-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-6-[[(2R,3R,4R,5S,6R)-3-acetamido-4-hydroxy-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-4,5-dihydroxyoxan-2-yl]oxy-6-[[(2S,3S,4S,5S,6R)-6-[[(2R,3R,4R,5S,6R)-3-acetamido-4-hydroxy-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-3,4,5-trihydroxyoxan-2-yl]oxymethyl]-3,5-dihydroxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-[[(2R,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxymethyl]oxan-3-yl]acetamide '
               'REASON: CORRECT Molecule has a large number of atoms (147) or '
               'repeating substructures (12)\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@@H]3[C@@H](OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]5CO)O[C@H]6[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]6CO[C@@H]7O[C@H]([C@@H](O)[C@@H](O)[C@@H]7O)C)O)[C@@H](O)[C@@H](O[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O[C@]%11(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%11)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]%10O)CO)[C@H](O)[C@H]9NC(=O)C)CO)CO)[C@@H]4O)O[C@@H]([C@@H](O)[C@@H]3O)CO)[C@@H]%12O[C@@H]([C@H](O)[C@H](O[C@]%13(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%13)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]%12O)CO '
               'NAME: '
               '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4S,5S,6R)-2-[(2R,3R,4S,5S,6S)-2-[[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-5-[(2S,3R,4S,5S,6R)-4-[(2S,4S,5R,6R)-5-acetamido-2-carboxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxan-2-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-[[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxymethyl]oxan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,5-dihydroxyoxan-4-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
               'acid REASON: CORRECT Molecule has a large number of atoms '
               '(172) or repeating substructures (13)\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'O1[C@H](O[C@H]2OC([C@@H](O)C(O)[C@H]2O)CO)C(O)C(O)[C@H](OC(=O)CCCCCCCCCCCCCCCCC)C1COC(=O)CCCCCCCCCCCCCCCCC '
               "NAME: alpha,alpha'-Trehalose 4,6-dioleate REASON: MISSED "
               'Molecule does not meet criteria for a macromolecule\n'
               ' * SMILES: '
               '[H+].[Br-].NCCCC[C@H](N)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(O)=O '
               'NAME: hexa-L-lysine monohydrobromide REASON: MISSED Molecule '
               'does not meet criteria for a macromolecule\n'
               ' * SMILES: '
               'O=C1NC=2C(=O)C(=CC(C2)=O)C(OC)C(CC(OC)C(OC)C(C=C(C(C(C=CC=C1C)C)OC(=O)N)C)C)C '
               'NAME: Macbecin 1 REASON: MISSED Molecule does not meet '
               'criteria for a macromolecule\n'
               ' * SMILES: '
               'O=C(O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@H](N)[C@H](O)C)[C@H](O)C1=CC=CC=C1)CC(C)C)C(C)C)[C@H]2NC(=NCC2)N)[C@H](O)C '
               'NAME: Iso-faulknamycin REASON: MISSED Molecule does not meet '
               'criteria for a macromolecule\n'
               ' * SMILES: '
               'O1[C@@](O[C@H]2O[C@@H]([C@@H](OC(=O)CCCCCCCCCCC)[C@H](O)[C@H]2O)CO)([C@@H](O)[C@H](O)[C@H]1CO)CO '
               'NAME: beta-D-Fructofuranosyl '
               '4-O-dodecanoyl-alpha-D-glucopyranoside REASON: MISSED Molecule '
               'does not meet criteria for a macromolecule\n'
               ' * SMILES: '
               'O=C(O[C@H]1[C@H](O[C@@H]2O[C@@H]([C@H](O)[C@@H]([C@H]2O)O)COC(=O)/C=C/C=C\\CCCCC)[C@H](OC([C@@H]1O)C3=C(O)C=C(O)C=C3CO)CO)/C(=C/C=C/CC(O)C(CC)C)/C '
               'NAME: Corynecandin REASON: MISSED Molecule does not meet '
               'criteria for a macromolecule\n'
               ' * SMILES: '
               'S1[C@H](CC(=O)N[C@H](C=2SC[C@H](C(N[C@@H](C1)CO)=O)N2)CC3=CC=CC=C3)CC(C)C '
               'NAME: Thioamycolamide C REASON: MISSED Molecule does not meet '
               'criteria for a macromolecule\n'
               ' * SMILES: '
               'S1C2=NC(=C1)C(=O)N[C@@H](C=3SC=C(N3)C(=O)N[C@H](C4=NC(C(N[C@H]2[C@H](CC)C)=O)=C(O4)C)C)CC(=O)O '
               'NAME: Microcyclamide GL546A REASON: MISSED Molecule does not '
               'meet criteria for a macromolecule\n'
               ' * SMILES: '
               'S1C/2=NC(=C1)C(=O)NC(C=3SC=C(N3)C(=O)NC(C=4SC=C(C5=C(C6=NC(C(NC(C(N\\C2=C/C)=O)C(O)C)=O)=CS6)C=CC(=N5)C=7SC=C(N7)C=8SC=C(N8)C(=O)N/C(/C(=O)NCC(=O)C)=C/C)N4)C(OC)C)C(O)(C)C '
               'NAME: YM-266184 REASON: MISSED Molecule does not meet criteria '
               'for a macromolecule\n'
               ' * SMILES: '
               'ClC1=C(OC)C=C2CC(=CC=C[C@@H](OC)[C@]3(OC)NC(O[C@H]([C@H]([C@H]4[C@]([C@H](CC(N(C1=C2)C)=O)OC(=O)CC)(O4)C)C)C3)=O)C '
               'NAME: Maytansinol-9-methyl ether-3-propionate REASON: MISSED '
               'Molecule does not meet criteria for a macromolecule\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'O1[C@H](O[C@H]2OC([C@@H](O)C(O)[C@H]2O)CO)C(O)C(O)[C@H](OC(=O)CCCCCCCCCCCCCCCCC)C1COC(=O)CCCCCCCCCCCCCCCCC',
                                      'name': "alpha,alpha'-Trehalose "
                                              '4,6-dioleate',
                                      'reason': 'Molecule does not meet '
                                                'criteria for a macromolecule'},
                                  {   'smiles': 'O=C1NC=2C(=O)C(=CC(C2)=O)C(OC)C(CC(OC)C(OC)C(C=C(C(C(C=CC=C1C)C)OC(=O)N)C)C)C',
                                      'name': 'Macbecin 1',
                                      'reason': 'Molecule does not meet '
                                                'criteria for a macromolecule'},
                                  {   'smiles': 'O1[C@@](O[C@H]2O[C@@H]([C@@H](OC(=O)CCCCCCCCCCC)[C@H](O)[C@H]2O)CO)([C@@H](O)[C@H](O)[C@H]1CO)CO',
                                      'name': 'beta-D-Fructofuranosyl '
                                              '4-O-dodecanoyl-alpha-D-glucopyranoside',
                                      'reason': 'Molecule does not meet '
                                                'criteria for a macromolecule'},
                                  {   'smiles': 'O=C(O[C@H]1[C@H](O[C@@H]2O[C@@H]([C@H](O)[C@@H]([C@H]2O)O)COC(=O)/C=C/C=C\\CCCCC)[C@H](OC([C@@H]1O)C3=C(O)C=C(O)C=C3CO)CO)/C(=C/C=C/CC(O)C(CC)C)/C',
                                      'name': 'Corynecandin',
                                      'reason': 'Molecule does not meet '
                                                'criteria for a macromolecule'},
                                  {   'smiles': 'S1[C@H](CC(=O)N[C@H](C=2SC[C@H](C(N[C@@H](C1)CO)=O)N2)CC3=CC=CC=C3)CC(C)C',
                                      'name': 'Thioamycolamide C',
                                      'reason': 'Molecule does not meet '
                                                'criteria for a macromolecule'},
                                  {   'smiles': 'S1C2=NC(=C1)C(=O)N[C@@H](C=3SC=C(N3)C(=O)N[C@H](C4=NC(C(N[C@H]2[C@H](CC)C)=O)=C(O4)C)C)CC(=O)O',
                                      'name': 'Microcyclamide GL546A',
                                      'reason': 'Molecule does not meet '
                                                'criteria for a macromolecule'},
                                  {   'smiles': 'S1C/2=NC(=C1)C(=O)NC(C=3SC=C(N3)C(=O)NC(C=4SC=C(C5=C(C6=NC(C(NC(C(N\\C2=C/C)=O)C(O)C)=O)=CS6)C=CC(=N5)C=7SC=C(N7)C=8SC=C(N8)C(=O)N/C(/C(=O)NCC(=O)C)=C/C)N4)C(OC)C)C(O)(C)C',
                                      'name': 'YM-266184',
                                      'reason': 'Molecule does not meet '
                                                'criteria for a macromolecule'},
                                  {   'smiles': 'ClC1=C(OC)C=C2CC(=CC=C[C@@H](OC)[C@]3(OC)NC(O[C@H]([C@H]([C@H]4[C@]([C@H](CC(N(C1=C2)C)=O)OC(=O)CC)(O4)C)C)C3)=O)C',
                                      'name': 'Maytansinol-9-methyl '
                                              'ether-3-propionate',
                                      'reason': 'Molecule does not meet '
                                                'criteria for a macromolecule'},
                                  {   'smiles': 'O1C(OCCCC(O)=O)C(O)C(O)C(O)C1C(O)=O',
                                      'name': '6-(3-carboxypropoxy)-3,4,5-trihydroxyoxane-2-carboxylic '
                                              'acid',
                                      'reason': 'Molecule does not meet '
                                                'criteria for a macromolecule'},
                                  {   'smiles': 'O(C1(OC(C(O)C1O)CO)CO)C2OC(C(O)C(O)C2O)COC3OC(C(O)C(O)C3O)COC(=O)/C(/C)=C/C',
                                      'name': '[6-[[6-[3,4-dihydroxy-2,5-bis(hydroxymethyl)oxolan-2-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-3,4,5-trihydroxyoxan-2-yl]methyl '
                                              '(E)-2-methylbut-2-enoate',
                                      'reason': 'Molecule does not meet '
                                                'criteria for a '
                                                'macromolecule'}],
    'attempt': 2,
    'reasoning': 'The previous program attempted to classify chemical entities '
                 'as macromolecules based on their SMILES strings by checking '
                 'for high molecular weight, large atom count, or repeating '
                 'substructures. While this criteria can often be associated '
                 'with macromolecules, the analysis failed to capture some '
                 'critical nuances; for example, some false negatives were due '
                 'to macromolecular features that the existing logic did not '
                 'capture, such as specific types of molecular repetition.\n'
                 '\n'
                 'To improve classification:\n'
                 '\n'
                 '1. **Consider Diverse Repeating Units:** Macromolecules can '
                 'include diverse repeating units, such as peptides, nucleic '
                 'acids, and polysaccharides. Improve substructural checks to '
                 'detect various types of repetitive units, which were not '
                 'adequately captured before.\n'
                 '\n'
                 '2. **Optimize Molecular Weight and Atom Count:** Although '
                 'high molecular weight is a hallmark of macromolecules, '
                 'improving detection might involve identifying known '
                 'macromolecular motifs or functional groups, which could be '
                 'small but significant.\n'
                 '\n'
                 '3. **Enhance Pattern Detection:** Employ SMARTS (SMiles '
                 'ARbitrary Target Specification) patterns to capture a '
                 'greater variety of repetitive substructures typical in '
                 'macromolecules.\n'
                 '\n'
                 '4. **False Negatives Analysis:** Check why certain known '
                 'macromolecules were mislabeled and adjust the pattern rules '
                 'or thresholds to account for these discrepancies.\n'
                 '\n'
                 'With these enhancements in mind, here is a revised version '
                 'of the program:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 721,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.8781973203410475,
    'f1': 0.9351491569390401,
    'accuracy': 0.8781973203410475,
    'negative_predictive_value': 0.0}