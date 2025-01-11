"""
Classifies: CHEBI:37739 glycerophospholipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    A glycerophospholipid must have a glycerol backbone with two fatty acid chains
    attached via ester linkages and a phosphate group ester-linked to a terminal carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophospholipid, False otherwise
        str: Reason for classification or failure
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone (O-C-C-C-O)
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)COP")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with phosphate found"

    # Look for two ester or ether groups (-O-C(=O)- or -O-C-)
    ester_or_ether_pattern = Chem.MolFromSmarts("[OX2][CX3](=O)[CX4]")
    ether_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    ester_matches = mol.GetSubstructMatches(ester_or_ether_pattern)
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    
    # Ensure there are at least two ester or ether groups
    if len(ester_matches) + len(ether_matches) < 2:
        return False, f"Found {len(ester_matches) + len(ether_matches)} ester or ether groups, need at least 2"

    # Check for phosphate group at terminal position
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No terminal phosphate group found"

    return True, "Contains glycerol backbone with two fatty acid chains and terminal phosphate group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37739',
                          'name': 'glycerophospholipid',
                          'definition': 'Any glycerolipid having a phosphate '
                                        'group ester-linked to a terminal '
                                        'carbon of the glycerol backbone.',
                          'parents': ['CHEBI:16247', 'CHEBI:35741'],
                          'xrefs': ['PMID:17393491'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O=C1C2=C(OC(=C1)C)C3=C(OC)C=C(OC)C=C3C(=C2O)C4=C5OC(=CC(C5=C(O)C=6C4=CC(OC)=CC6OC)=O)C',
                                     'name': 'Isonigerone',
                                     'reason': 'No glycerol backbone with '
                                               'phosphate found'},
                                 {   'smiles': 'C(=C\\C/C=C\\CCCC(NC)=O)\\C/C=C\\C/C=C\\CCCCC',
                                     'name': 'N-methyl arachidonoyl amine',
                                     'reason': 'No glycerol backbone with '
                                               'phosphate found'},
                                 {   'smiles': 'OC1(O)[C@]23N(CC1)C(N(O)[C@H]([C@@]3(N=C(N2)N)[H])CO)=N',
                                     'name': 'Decarbamoylneosaxitoxin',
                                     'reason': 'No glycerol backbone with '
                                               'phosphate found'},
                                 {   'smiles': 'S([C@H]1N(C(=O)/C(=C\\C2=CC=C(OCC=C(C)C)C=C2)/N(C1=O)C)C)C',
                                     'name': 'Fusaperazine F',
                                     'reason': 'No glycerol backbone with '
                                               'phosphate found'},
                                 {   'smiles': 'Oc1ccccc1I',
                                     'name': '2-iodophenol',
                                     'reason': 'No glycerol backbone with '
                                               'phosphate found'},
                                 {   'smiles': 'O=C1N(CC(=O)N[C@H](C(=O)O[C@H]([C@@H](C(NC(C=C1)=C)=O)C)C(CCCCCCCCCCCCCC)C)C(O)C(=O)N)C',
                                     'name': 'Rakicidin H',
                                     'reason': 'No glycerol backbone with '
                                               'phosphate found'},
                                 {   'smiles': 'O(C1=C(OC)C=C(C2=C(O)C(OC)=C(C3=CC=CC=C3)C=C2OC)C=C1)C/C=C(/CO)\\C',
                                     'name': 'Prenylterphenyllin F',
                                     'reason': 'No glycerol backbone with '
                                               'phosphate found'},
                                 {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)[C@@H](O)[C@@H]1OC[C@H]3O[C@@H](OC[C@H]4O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]4O)[C@H](O)[C@@H](O)[C@@H]3O)CO[C@@H]5O[C@@H]([C@@H](O)[C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)[C@H]5O)CO[C@@H]7O[C@@H]([C@@H](O)[C@H](O)[C@H]7O)CO',
                                     'name': '(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5R,6R)-6-[[(2R,3R,4S,5R,6R)-3,5-Dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-2-yl]oxymethyl]-3,5-dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-3,4,5-trihydroxyoxan-2-yl]oxymethyl]oxane-2,3,4,5-tetrol',
                                     'reason': 'No glycerol backbone with '
                                               'phosphate found'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)O)[C@@H]4O[C@@H]([C@@H](O)[C@H](O[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O)CO)[C@H](O)[C@@H]6NC(=O)C)CO)CO)[C@@H]4O)CO[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO)[C@H](O)[C@H]9NC(=O)C)CO)CO',
                                     'name': 'Gal2GlcNAc2Man3GlcNAcFucGlcNAc',
                                     'reason': 'No glycerol backbone with '
                                               'phosphate found'},
                                 {   'smiles': 'CN1CC2(CCN(CC2)S(=O)(=O)C3=CN(C=N3)C)C4=C([C@@H]1CO)NC5=C4C=CC(=C5)OC',
                                     'name': "[(1R)-7-methoxy-2-methyl-1'-[(1-methyl-4-imidazolyl)sulfonyl]-1-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]yl]methanol",
                                     'reason': 'No glycerol backbone with '
                                               'phosphate found'}],
    'sample_false_negatives': [   {   'smiles': 'O[C@@H]1C/C(=C\\C=C\\2/[C@]3([C@@]([C@](CC3)([C@H](C)C=C[C@@H](C(C)C)C)[H])(CCC2)C)[H])/C([C@@H](O)C1)=C',
                                      'name': '1-octadecyl lysophosphatidic '
                                              'acid',
                                      'reason': 'No glycerol backbone with '
                                                'phosphate found'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCSC[C@H](COP(O)(=O)OCC[N+](C)(C)C)OP(O)(=O)CCCCCCCCCCCCCCCC',
                                      'name': '1-S-hexadecyl-2-O-[hexadecyl(hydroxy)phosphoryl]-1-thio-sn-glycero-3-phosphocholine',
                                      'reason': 'No glycerol backbone with '
                                                'phosphate found'},
                                  {   'smiles': '[C@](COC(=O)CCCCCCCCCCCCCCC)(COP([O-])(=O)CC[N+](C)(C)C)(OC(=O)CCCCCCC/C=C\\CCCCCCCC)[H]',
                                      'name': 'PnC(16:0/18:1(9Z))',
                                      'reason': 'No terminal phosphate group '
                                                'found'},
                                  {   'smiles': 'C(CCCCCCCCCCCC)CCCSC[C@@H](NC(CCCCCCCCCCCCCCC)=O)COP(OCC[N+](C)(C)C)(=O)[O-]',
                                      'name': 'thioetheramide PC',
                                      'reason': 'No glycerol backbone with '
                                                'phosphate found'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCC(=O)OCC(COP(O)(=O)CC[N+](C)(C)C)OC(=O)CCCCCCC\\C=C/CCCCCCCC',
                                      'name': '1-hexadecanoyl-2-[(Z)-octadec-9-enoyl]-sn-glycerol '
                                              '3-(2-trimethylaminoethylphosphonate)',
                                      'reason': 'No terminal phosphate group '
                                                'found'},
                                  {   'smiles': 'O1C(C1)COC(=O)CCCC=CCC=CCC=CCC=CCCCCC',
                                      'name': 'O-Arachidonoyl Glycidol',
                                      'reason': 'No glycerol backbone with '
                                                'phosphate found'},
                                  {   'smiles': 'C(CCCCCCC/C=C\\CCCCCCCC)(OC[C@@H](COP(O)(=S)O)OC)=O',
                                      'name': '(2S)-1-oleoyl-2-methylglycero-3-phosphothionate',
                                      'reason': 'No terminal phosphate group '
                                                'found'},
                                  {   'smiles': 'P(OC(COC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCC=CCCCCCCCC)(OCC[NH3+])([O-])=O',
                                      'name': '1-palmitoyl-3-oleoyl-sn-glycero-2-phosphoethanolamine',
                                      'reason': 'No glycerol backbone with '
                                                'phosphate found'}],
    'attempt': 0,
    'reasoning': 'To classify a molecule as a glycerophospholipid using its '
                 'SMILES representation, we need to identify specific '
                 'structural features. A glycerophospholipid consists of a '
                 'glycerol backbone, to which two fatty acid chains are '
                 'attached via ester bonds, and a phosphate group ester-linked '
                 'to a terminal carbon of the glycerol.\n'
                 '\n'
                 "Here's how we can develop a step-by-step approach to this "
                 'classification:\n'
                 '\n'
                 '1. **Validate the SMILES string**: First, we need to check '
                 'if the SMILES string is valid and can be converted into a '
                 'molecular structure.\n'
                 '\n'
                 '2. **Identify the Glycerol Backbone**: Using a substructural '
                 'match, we should detect a three-carbon chain with two ester '
                 'or ether linkages at the C1 and C2 positions and a hydroxyl '
                 'or similar group at C3 to which the phosphate can be '
                 'linked. \n'
                 '\n'
                 '3. **Detect the Phosphate Group**: A glycerophospholipid '
                 'must have a phosphate group attached to the terminal carbon '
                 'of the glycerol backbone. The SMARTS pattern to find a '
                 'phosphate group -linked to oxygen is necessary.\n'
                 '\n'
                 '4. **Verify the Structure**: Ensure that there are fatty '
                 'acid chains linked to the glycerol (possibly as ester or '
                 'ether) and the phosphate group is present as described.\n'
                 '\n'
                 'Below is a Python code using RDKit to implement this logic:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 3867,
    'num_false_positives': 100,
    'num_true_negatives': 2879,
    'num_false_negatives': 8,
    'num_negatives': None,
    'precision': 0.9747920342828333,
    'recall': 0.9979354838709678,
    'f1': 0.9862280030604438,
    'accuracy': 0.9842427779398891,
    'negative_predictive_value': 0.99722895739522}