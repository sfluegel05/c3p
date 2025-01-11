"""
Classifies: CHEBI:28892 ganglioside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    A ganglioside is a glycosphingolipid containing one or more sialic acids 
    linked to an oligosaccharide chain attached to a ceramide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight - gangliosides are large molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800:
        return False, "Molecular weight too low for ganglioside"

    # Look for ceramide backbone
    # Sphingosine pattern (long chain base with OH and NH)
    sphingosine_pattern = Chem.MolFromSmarts("[CH2]-[CH]-[CH](-[OH])-[CH](-[NH]-[C](=O))-[CH2]-O")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No ceramide backbone found"

    # Look for sialic acid (Neu5Ac) pattern
    # Core structure of sialic acid with carboxyl group
    sialic_pattern = Chem.MolFromSmarts("[C](=O)[O;H,-]-[C]1-[C]-[C](-[NH]-[C](=O)-[CH3])-[C]-[C](-[OH])-[C]1")
    sialic_matches = mol.GetSubstructMatches(sialic_pattern)
    if len(sialic_matches) < 1:
        return False, "No sialic acid (Neu5Ac) found"

    # Look for oligosaccharide pattern
    # Basic sugar rings with OH groups
    sugar_pattern = Chem.MolFromSmarts("[C]1-[O]-[C]-[C]-[C]-[C]-1")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) < 3:
        return False, "Insufficient oligosaccharide chain"

    # Check for glycosidic linkages
    glycosidic_pattern = Chem.MolFromSmarts("[C]-[O]-[C]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if len(glycosidic_matches) < 4:
        return False, "Insufficient glycosidic linkages"

    # Count key atoms to verify overall composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 30:
        return False, "Too few carbons for ganglioside structure"
    if o_count < 15:
        return False, "Too few oxygens for ganglioside structure"
    if n_count < 2:
        return False, "Too few nitrogens for ganglioside structure"

    # Verify presence of long chain fatty acid
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "Missing long chain fatty acid"

    return True, "Contains ceramide backbone, sialic acid(s), and oligosaccharide chain with correct linkages"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:28892',
                          'name': 'ganglioside',
                          'definition': 'A molecule composed of a '
                                        'glycosphingolipid (ceramide and '
                                        'oligosaccharide) with one or more '
                                        'sialic acids linked on the sugar '
                                        'chain.',
                          'parents': ['CHEBI:17761', 'CHEBI:36526'],
                          'xrefs': [   'MetaCyc:Gangliosides',
                                       'PMID:16158191',
                                       'PMID:2088646',
                                       'PMID:38623278',
                                       'PMID:38848944',
                                       'PMID:38887845',
                                       'PMID:39092231',
                                       'Wikipedia:Ganglioside'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'Cc1cc(C)c2Cc3c(C)cc(C)c(Cc4c(C)cc(C)c(Cc5c(C)cc(C)c(Cc1c2C)c5C)c4C)c3C',
                                     'name': '4,6,10,12,16,18,22,24,25,26,27,28-dodecamethylcalix[4]arene',
                                     'reason': 'Molecular weight too low for '
                                               'ganglioside'},
                                 {   'smiles': 'O=C1C(=C2C=C3[C@]([C@@H](C(C)C)[C@@H]([C@H]3O)OC(=O)C)(C)CC[C@]2(C)CC1)COC(=O)C',
                                     'name': 'Dahliane E',
                                     'reason': 'Molecular weight too low for '
                                               'ganglioside'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C',
                                     'name': 'Glu-Thr-Leu',
                                     'reason': 'Molecular weight too low for '
                                               'ganglioside'},
                                 {   'smiles': 'CCOc1ccc(NC(=O)C(C)O)cc1',
                                     'name': 'p-Lactophenetide',
                                     'reason': 'Molecular weight too low for '
                                               'ganglioside'},
                                 {   'smiles': 'O=C1NCC=2C1=C3C(N([C@@H]4O[C@H]([C@@H](O)[C@H]([C@H]4OC)O)C)C5=C3C=CC=C5)=C6NC7=C(C26)C=CC=C7',
                                     'name': "3'-epi-5'-methoxy-K252d",
                                     'reason': 'Molecular weight too low for '
                                               'ganglioside'},
                                 {   'smiles': 'O=C(OC1C(O)C(OC(C1O)CO)OCC2OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(O)C(C2O)O)CCCCCCCCCCCCCCC',
                                     'name': '[(2S)-2-hexadecanoyloxy-3-[(2S,3R,4S,5S,6R)-6-[[(2S,3R,4S,5S,6R)-4-hexadecanoyloxy-3,5-dihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl]oxymethyl]-3,4,5-trihydroxy-tetrahydropyran-2-yl]oxy-propyl] '
                                             '(9E,12E,15E)-octadeca-9,12,15-trienoate',
                                     'reason': 'No ceramide backbone found'},
                                 {   'smiles': 'O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@@H](O)[C@]([C@@H]3O)(O)C)C',
                                     'name': 'Altersolanol G',
                                     'reason': 'Molecular weight too low for '
                                               'ganglioside'},
                                 {   'smiles': '[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(=O)CO)C([O-])=O)[C@H](O)[C@H](O)CO',
                                     'name': 'N-glycoloyl-alpha-neuraminate',
                                     'reason': 'Molecular weight too low for '
                                               'ganglioside'},
                                 {   'smiles': 'OC(C(O)C/C=C\\C/C=C\\C/C=C\\CC)C/C=C\\C/C=C\\CCC(O)=O',
                                     'name': '10,11-DiHDPE',
                                     'reason': 'Molecular weight too low for '
                                               'ganglioside'},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'Molecular weight too low for '
                                               'ganglioside'}],
    'sample_false_negatives': [   {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O[C@@H]3O[C@H](CO)[C@H](O)[C@H](O[C@@H]4O[C@H](CO)[C@H](O)[C@H](O[C@@]5(C[C@H](O)[C@@H](NC(C)=O)[C@@H](O5)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]4O)[C@H]3NC(C)=O)[C@H](O[C@@]3(C[C@H](O)[C@@H](NC(C)=O)[C@@H](O3)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]2O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC',
                                      'name': "alpha-Neu5Ac-(2->3)-beta-D-Gal-(1->3)-beta-D-GalNAc-(1->4)-[alpha-Neu5Ac-(2->3)]-beta-D-Gal-(1->4)-beta-D-Glc-(1<->1')-Cer(d18:1/22:0)",
                                      'reason': 'No ceramide backbone found'},
                                  {   'smiles': '[C@H]1(O[C@@H]([C@@H]([C@@H]([C@H]1O)O[C@@]23O[C@]([C@@H]([C@H](C2)O)NC3=O)([C@@H]([C@@H](CO)O)O)[H])O)CO)O[C@H]4[C@@H]([C@H]([C@@H](O[C@@H]4COS(=O)(=O)O)O[C@@H]5[C@H]([C@@H](O[C@@H]([C@@H]5O)CO)O[C@H]6[C@@H]([C@H]([C@@H](O[C@@H]6CO)OC[C@@H]([C@@H](/C=C/CCCCCCCCCCCCC)O)NC(=O)CCCCCCCCCCCCCCCCC)O)O)O)NC(=O)C)O',
                                      'name': "alpha-Neu5d5N;1,5-lactam-(2->3)-beta-D-Gal-(1->4)-beta-D-GlcNAc6S-(1->3)-beta-D-Gal-(1->4)-beta-D-Glc-(1<->1')-Cer(d18:1/18:0)",
                                      'reason': 'No ceramide backbone found'},
                                  {   'smiles': '[C@@]1(O[C@]([C@H](NC(=O)C)[C@H](C1)O)([C@@H]([C@H](O)CO)O)[H])(C(O)=O)O[C@H]2[C@H]([C@H](O[C@H]([C@@H]2O)O[C@@H]3[C@H](O[C@@H](OC[C@@H]([C@@H](/C=C/CCCCCCCCCCCCC)O)NC(=O)CCCCCCCCCCCCCCC)[C@@H]([C@H]3O)O)CO)CO)O[C@H]4[C@@H]([C@H]([C@@H](O)[C@H](O4)CO)O)NC(C)=O',
                                      'name': 'ganglioside GM2 (16:0)',
                                      'reason': 'No ceramide backbone found'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O[C@@]3(C[C@H](O)[C@@H](NC(C)=O)[C@@H](O3)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]2O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC',
                                      'name': "alpha-Neu5Ac-(2->3)-beta-D-Gal-(1->4)-beta-D-Glc-(1<->1')-Cer(d18:1/18:0)",
                                      'reason': 'No ceramide backbone found'},
                                  {   'smiles': 'O([C@@H]1[C@H]([C@@H](O[C@@H]([C@@H]1O)CO)O[C@H]2[C@@H]([C@H]([C@@H](O[C@@H]2CO)OC[C@@H]([C@@H](/C=C/CCCCCCCCCCCCC)O)NC(CCCCCCCCCCCCCCCCC)=O)O)O)O)[C@]3(O[C@]([C@@H]([C@H](C3)O)NC(CO)=O)([C@@H]([C@H](O)CO)O)[H])C(=O)O',
                                      'name': 'N-glycolylneuraminyllactosylceramide',
                                      'reason': 'No ceramide backbone found'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O[C@@H]3O[C@H](CO)[C@H](O)[C@H](O[C@@H]4O[C@H](CO)[C@H](O)[C@H](O[C@@]5(C[C@H](O)[C@@H](NC(C)=O)[C@@H](O5)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]4O)[C@H]3NC(C)=O)[C@H](O[C@@]3(C[C@H](O)[C@@H](NC(C)=O)[C@@H](O3)[C@H](O[C@@]3(C[C@H](O)[C@@H](NC(C)=O)[C@@H](O3)[C@H](O[C@@]3(C[C@H](O)[C@@H](NC(C)=O)[C@@H](O3)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)CO)C(O)=O)[C@H](O)CO)C(O)=O)[C@H]2O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC',
                                      'name': "alpha-Neu5Ac-(2->3)-beta-D-Gal-(1->3)-beta-D-GalNAc-(1->4)-[alpha-Neu5Ac-(2->8)-alpha-Neu5Ac-(2->8)-alpha-Neu5Ac-(2->3)]-beta-D-Gal-(1->4)-beta-D-Glc-(1<->1')-Cer(d18:1/22:0)",
                                      'reason': 'No ceramide backbone found'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@H](O[C@@]3(C[C@H](O)[C@@H](NC(C)=O)[C@@H](O3)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]2O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC',
                                      'name': "alpha-Neu5Ac-(2->3)-beta-D-Gal-(1->4)-beta-D-Glc-(1<->1')-Cer(d18:1/23:0)",
                                      'reason': 'No ceramide backbone found'},
                                  {   'smiles': '[C@@]1(O[C@]([C@@H]([C@H](C1)O)NC(C)=O)([C@@H]([C@@H](CO)O)O)[H])(C(O)=O)O[C@@H]([C@H]([C@@]2(O[C@](O[C@H]3[C@H]([C@H](O[C@H]([C@@H]3O)O[C@@H]4[C@H](O[C@H]([C@@H]([C@H]4O)O)OC[C@@H]([C@@H](/C=C/CCCCCCCCCCCCC)O)NC(CCCCCCCCCCCCCCCCC)=O)CO)CO)O[C@H]5[C@@H]([C@H]([C@H]([C@H](O5)CO)O)O[C@@H]6O[C@@H]([C@@H]([C@@H]([C@H]6O)O)O)CO)NC(C)=O)(C[C@H](O)[C@H]2NC(=O)C)C(O)=O)[H])O)CO',
                                      'name': "beta-Gal-(1->3)-beta-GalNAc-(1->4)-[alpha-Neu5Ac-(2->8)-alpha-Neu5Ac-(2->3)]-beta-Gal-(1->4)-beta-Glc-(1<->1')-Cer(d18:1/18:0)",
                                      'reason': 'No ceramide backbone found'},
                                  {   'smiles': '[C@@H]1([C@@H]([C@H]([C@H]([C@H](O1)CO)O[C@H]2[C@@H]([C@H]([C@H]([C@H](O2)CO)O)O[C@@H]3O[C@@H]([C@@H]([C@@H]([C@H]3O)O[C@]4(O[C@]([C@@H]([C@H](C4)O)NC(C)=O)([C@@H]([C@@H](CO)O)O)[H])C(O)=O)O)CO)NC(=O)C)O[C@]5(O[C@]([C@H](NC(=O)C)[C@H](C5)O)([C@@H]([C@H](O)CO)O)[H])C(=O)O)O)O[C@H]6[C@@H]([C@H]([C@@H](O[C@@H]6CO)OC[C@@H]([C@@H](/C=C/CCCCCCCCCCCCC)O)NC(=O)CCCCCCCCCCCCCCC)O)O',
                                      'name': "alpha-Neu5Ac-(2->3)-beta-D-Gal-(1->3)-beta-D-GalNAc-(1->4)-[alpha-Neu5Ac-(2->3)]-beta-D-Gal-(1->4)-beta-D-Glc-(1<->1')-Cer(d18:1/16:0)",
                                      'reason': 'No ceramide backbone found'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O[C@@H]3O[C@H](CO)[C@H](O)[C@H](O[C@@H]4O[C@H](CO)[C@H](O)[C@H](O)[C@H]4O)[C@H]3NC(C)=O)[C@H](O[C@@]3(C[C@H](O)[C@@H](NC(C)=O)[C@@H](O3)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]2O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC',
                                      'name': 'ganglioside GM1',
                                      'reason': 'No ceramide backbone found'}],
    'attempt': 0,
    'reasoning': "I'll help create a program to classify gangliosides. "
                 'Gangliosides are complex glycosphingolipids containing one '
                 'or more sialic acids (Neu5Ac) linked to an oligosaccharide '
                 "chain attached to a ceramide. Let's break down the key "
                 'structural features:\n'
                 '\n'
                 '1. Must contain ceramide (sphingosine + fatty acid)\n'
                 '2. Must have one or more sialic acids (Neu5Ac)\n'
                 '3. Must have oligosaccharide chain (usually containing Gal, '
                 'GalNAc, Glc)\n'
                 '4. Sialic acids must be linked to the sugar chain\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 142250,
    'num_false_negatives': 50,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.999648629655657,
    'negative_predictive_value': 0.999648629655657}