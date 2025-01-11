"""
Classifies: CHEBI:15841 polypeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide based on its SMILES string.
    Polypeptides are characterized by sequences of peptide bonds (-CO-NH-) and amino acid structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # The peptide bond pattern may need refinement to capture all necessary forms, including cyclized forms
    peptide_bond_pattern = Chem.MolFromSmarts("N[C;R0][C;R0](=O)")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) < 5:  # Assuming at least 5 peptide bonds indicate a polypeptide
        return False, f"Detected {len(peptide_bond_matches)} peptide bonds, need at least 5 for a polypeptide"

    # Count the number of nitrogen atoms generally found in polypeptides
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if nitrogen_count < 10:  # Expecting at least 10 nitrogen atoms typically found in short to long peptides
        return False, f"Detected {nitrogen_count} nitrogen atoms, too low for a polypeptide"

    # Check molecular weight to ensure we capture larger polypeptide molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:  # Molecular weight of 500 Dalton is typically a minimum for a polypeptide
        return False, f"Molecular weight too low ({mol_wt} Da) for a polypeptide"

    # Make sure there's a continuous backbone structure, i.e., no isolated segments of peptide-like groups
    molecular_fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(molecular_fragments) > 1:
        return False, "Polypeptide structure should not be fragmented"

    return True, "Contains characteristic peptide bonds and nitrogen-rich backbone indicative of a polypeptide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:15841',
                          'name': 'polypeptide',
                          'definition': 'A peptide containing ten or more '
                                        'amino acid residues.',
                          'parents': ['CHEBI:16670', 'CHEBI:33839'],
                          'xrefs': ['KEGG:C00403'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.05660377358490566 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'CC[C@H](C)[C@H](N)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](C)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](C)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H]([C@@H](C)CC)C(=O)N[C@@H](CC(C)C)C(N)=O '
               'NAME: mastoparan REASON: CORRECT Contains peptide bonds and '
               'polypeptide backbone structure with sufficient molecular '
               'complexity\n'
               ' * SMILES: '
               'O=C(NCC(=O)N[C@H](C(=O)N[C@H](CO)CC(C)C)[C@H](CC)C)C(NC(=O)[C@@H](NC(=O)CNC(=O)CNC(=O)C(NC(=O)[C@@H](NC(=O)CNC(=O)C(NC(=O)CCCCCCC)(C)C)CC(C)C)(C)C)CC(C)C)(C)C '
               'NAME: Trichogin A IV REASON: CORRECT Contains peptide bonds '
               'and polypeptide backbone structure with sufficient molecular '
               'complexity\n'
               ' * SMILES: '
               'O=C(O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CC(O)CCCCCCCCC)CC(C)C)CC(C)C)CCC(=O)N)CC(C)C)C(O)C)C(CC)C)CC(C)C)CC(C)C '
               'NAME: Syringafactin F REASON: CORRECT Contains peptide bonds '
               'and polypeptide backbone structure with sufficient molecular '
               'complexity\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'SCC(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(N)CCCCN)CC(=O)N)CC1=CC=C(O)C=C1)CC(=O)N)C(C)C)CC=2NC=NC2)C(=O)NC(C(=O)N)[C@H](O)C '
               'NAME: Leucocin C-TA33a REASON: MISSED Molecular complexity too '
               'low for a polypeptide\n'
               ' * SMILES: '
               'O=C1NCCC(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CC(O)CCC(CC)C)CCN)C(O)C)CCN)C(=O)NC(C(=O)NC(C(=O)NC(C(NC(C(NC(C(NC1C(O)C)=O)CCN)=O)CCN)=O)CC(C)C)CC2=CC=CC=C2)CCN '
               'NAME: Polymyxin B6 REASON: MISSED Molecular complexity too low '
               'for a polypeptide\n'
               ' * SMILES: '
               'SCC1NC(=O)C2NC(=O)C3NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(N)C(C)C)CCC(=O)N)C)C(OC(=O)CC(NC(=O)C(CCC(NCCCC3)=O)NC(CNC(C(NC(C(NC(C(NC(C(NC1=O)CC(=O)O)=O)C(C)C)=O)CCC(NCCCC2)=O)=O)C(O)C)=O)=O)C(=O)NC(C(=O)NC(C(=O)NCC(=O)NC(C(=O)O)CCC(=O)O)CCCCN)CCC(=O)N)C '
               'NAME: Marinomonasin REASON: MISSED Molecular complexity too '
               'low for a polypeptide\n'
               ' * SMILES: '
               'O=C(NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC(C)C)CCC(=O)N)(C)C)CC(C)C)CCC(=O)N)C)C)(C)C)CCC(=O)N)(C)C)CC(C)C)(C)C)C(NC(=O)C(NC(=O)C(NC(=O)C)CC=1C2=C(C=CC=C2)NC1)C)(C)C '
               'NAME: Ampullosporin E2 REASON: MISSED Molecular complexity too '
               'low for a polypeptide\n'
               ' * SMILES: '
               'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC2=CC=CC=C2)CCC(=O)N)CCC(=O)N)(C)C)(C)C)C(C)C)CCC1)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)(C)C)C)(C)C)C)(C)C)C)CCC(=O)N)(C)C)C(CC)C)(C)C)CC(C)C)(C)C '
               'NAME: Trichocellin-A-VII REASON: MISSED Molecular complexity '
               'too low for a polypeptide\n'
               ' * SMILES: '
               'S=C(N[C@H](C(=O)N[C@@H]1C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@@H](CC2=CC=C(OC)C=C2)C(N[C@H](C(NC=CS[C@@H]1C)=O)C(O)C=3N(C=[N+](C3)C)C)=O)C)C(C)C)C)[C@@H](NC(=S)[C@@H](NC(=O)[C@@H](NC(=S)[C@@H](NC(=S)[C@@H](NC(=O)C(=O)C)C(C)C)CCSC)C)C)C(CC)C '
               'NAME: Thiostreptamide S4 REASON: MISSED Molecular complexity '
               'too low for a polypeptide\n'
               ' * SMILES: '
               '[H]N[C@H]1CSSC[C@H](NC([C@H]([C@@H](C)O)NC([C@@H](NC([C@H](CC(C)C)NC([C@@H](NC([C@@H](NC1=O)CO)=O)CC(N)=O)=O)=O)CO)=O)=O)C(N[C@@H](CC(C)C)C(N[C@@H](CC(C)C)C(N[C@@H](CO)C(N[C@@H](C)C(N[C@@H](CC2=CC=C(O)C=C2)C(N[C@@H](CC3=CNC4=C3C=CC=C4)C(N[C@@H](CCCCN)C(N[C@@H](CC(O)=O)C(N[C@@H](CC(C)C)C(N[C@@H](CC(N)=O)C(N[C@@H](CC(N)=O)C(N[C@@H](CC5=CC=C(O)C=C5)C(N[C@@H](CC6=CNC=N6)C(N[C@@H](CCCNC(N)=N)C(N[C@@H](CC7=CC=CC=C7)C(N[C@@H](CO)C(NCC(N[C@@H](CCSC)C(NCC(N[C@@H](CC8=CC=CC=C8)C(NCC(N9CCC[C@H]9C(N[C@@H](CCC(O)=O)C(N[C@@H]([C@@H](C)O)C(N%10CCC[C@H]%10C(N)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O '
               'NAME: calcitonin (pork natural) REASON: MISSED Molecular '
               'complexity too low for a polypeptide\n'
               ' * SMILES: '
               'O=C(O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CCCCC)C(C)C)CC1=CC=CC=C1)CCC(=O)O)CC(C)C)C)CC(C)C '
               'NAME: Holrhizin P REASON: MISSED Molecular complexity too low '
               'for a polypeptide\n'
               ' * SMILES: '
               'CC[C@H](C)[C@@H](C(=O)N[C@@]1(CCCC=CCCC[C@](NC(=O)[C@@H](NC(=O)[C@@H]2CCCN2C(=O)[C@@H](NC1=O)CCSC)CC(C)C)(C)C(=O)N[C@@H](C)C(=O)N3CCC[C@H]3C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC4=CC=C(C=C4)O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N)C)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](C(C)C)NC(=O)CCNC(=O)C '
               'NAME: LSM-37192 REASON: MISSED Molecular complexity too low '
               'for a polypeptide\n'
               ' * SMILES: '
               'S1[C@H](NC(C1)C(=O)N[C@H](C(=O)O)CCSC)[C@H]2N(C(=O)[C@@H](NC(=O)/C(/NC(=O)[C@@H](NC(=O)/C(/NC(=O)[C@@H](NC(=O)/C(/NC(=O)N[C@H](C(=O)O)CC3=CC=CC=C3)=C/C)CO)=C/C)CO)=C/C)CCCN=C(N)N)CCC2 '
               'NAME: Chitinimide F REASON: MISSED Molecular complexity too '
               'low for a polypeptide\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'O=C1NCCC(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CC(O)CCC(CC)C)CCN)C(O)C)CCN)C(=O)NC(C(=O)NC(C(=O)NC(C(NC(C(NC(C(NC1C(O)C)=O)CCN)=O)CCN)=O)CC(C)C)CC2=CC=CC=C2)CCN',
                                      'name': 'Polymyxin B6',
                                      'reason': 'Detected 3 peptide bonds, '
                                                'need at least 5 for a '
                                                'polypeptide'},
                                  {   'smiles': 'S=C(N[C@H](C(=O)N[C@@H]1C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@@H](CC2=CC=C(OC)C=C2)C(N[C@H](C(NC=CS[C@@H]1C)=O)C(O)C=3N(C=[N+](C3)C)C)=O)C)C(C)C)C)[C@@H](NC(=S)[C@@H](NC(=O)[C@@H](NC(=S)[C@@H](NC(=S)[C@@H](NC(=O)C(=O)C)C(C)C)CCSC)C)C)C(CC)C',
                                      'name': 'Thiostreptamide S4',
                                      'reason': 'Detected 3 peptide bonds, '
                                                'need at least 5 for a '
                                                'polypeptide'},
                                  {   'smiles': 'O=C(O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CCCCC)C(C)C)CC1=CC=CC=C1)CCC(=O)O)CC(C)C)C)CC(C)C',
                                      'name': 'Holrhizin P',
                                      'reason': 'Detected 6 nitrogen atoms, '
                                                'too low for a polypeptide'},
                                  {   'smiles': 'S1CC(NC(=O)C(NC(=O)C2NC(=O)CNC(=O)C3CC(O)CN3C(C(C(SC2)C)NC(=O)C4NC(=O)C(NC(=O)C(NC(=O)C(NC(C(CSC4)NC(=O)/C(/NC(=O)C(N)C(CC)C)=C/C)=O)C(CC)C)=C)CC(C)C)=O)C(O)C)C(=O)NC(C(=O)NCC(=O)NCC(=O)NCC(NC5C(NC(C(NC(C1)C(=O)NC6C(=O)NC(C(=O)NC(CSC5)C(NC=CSC6)=O)CC7=CC=CC=C7)=O)CC(=O)N)=O)=O)C',
                                      'name': 'Daspyromycin B',
                                      'reason': 'Detected 3 peptide bonds, '
                                                'need at least 5 for a '
                                                'polypeptide'},
                                  {   'smiles': 'O=C1N[C@@H]2C(=O)N[C@H](C(=O)O)CC3=CC=C(OC=4C=C(C2)C=C(C5=C(OC)C=CC(C[C@@H](C(N[C@H]1[C@H](CC)C)=O)N)=C5)C4O)C=C3',
                                      'name': 'Cittilin A',
                                      'reason': 'Detected 0 peptide bonds, '
                                                'need at least 5 for a '
                                                'polypeptide'},
                                  {   'smiles': 'O=C(O)C(NC(=O)C(NC(=O)[C@@H](NC(=O)CC(O)CC(O)C(NC(=O)CCCCCCC/C=C\\CCCCCC)CC=1NC=NC1)C)CCCN(O)C(=O)C)CCCN(O)C(=O)C',
                                      'name': 'Pseudoalteropeptide A',
                                      'reason': 'Detected 3 peptide bonds, '
                                                'need at least 5 for a '
                                                'polypeptide'},
                                  {   'smiles': 'CC[C@H](C)[C@H](N)C1=N[C@@H](CS1)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](CCC(O)=O)C(=O)N[C@@H]([C@@H](C)CC)C(=O)N[C@H]1CCCCNC(=O)[C@H](CC(N)=O)NC(=O)[C@@H](CC(O)=O)NC(=O)[C@H](Cc2c[nH]cn2)NC(=O)[C@@H](Cc2ccccc2)NC(=O)[C@@H](NC(=O)[C@@H](CCCN)NC1=O)[C@@H](C)CC',
                                      'name': 'bacitracin A',
                                      'reason': 'Detected 3 peptide bonds, '
                                                'need at least 5 for a '
                                                'polypeptide'},
                                  {   'smiles': 'S1C(=NC(=C1)C(=O)NCC(=O)NC(C(=O)NCC(=O)NC2C[N+](C=3C=C(O)C(=CC3C2)O)(C)C)C)C(NC(=O)C(O)(C(=O)O)CCNC(=O)C(NC(=O)C(O)(C(=O)O)CC(O)CCCCCCC)C)C',
                                      'name': 'Fabrubactin B',
                                      'reason': 'Detected 4 peptide bonds, '
                                                'need at least 5 for a '
                                                'polypeptide'},
                                  {   'smiles': 'ClC1=C2OC3=C(O[C@H]4O[C@H]([C@H](O)[C@@H]([C@@H]4N)O)CO)C=5OC6=CC=C(CC7NC(=O)C(N)C=8C=CC(=C(OC=9C=C(C(C(NC(C(=C3)C5)C(=O)NC%10C(NC(C(C(=C1)C=C2)O[C@H]%11O[C@@H]([C@@H](O)[C@@H]([C@H]%11NC(=O)C)O)CO)C(=O)NC(C%12=C(C(O[C@H]%13OC([C@@H](O)[C@@H]([C@@H]%13O)O)CO)=CC(=C%12)O)C=%14C=C%10C=CC%14O)C(=O)O)=O)=O)NC7=O)C=C(O)C9)C8)O)C=C6',
                                      'name': 'Deacyl-GP1416',
                                      'reason': 'Detected 0 peptide bonds, '
                                                'need at least 5 for a '
                                                'polypeptide'},
                                  {   'smiles': 'S=C(NC1C(=O)N[C@@H]2C(=O)NC3C(=O)NC(C(=O)N[C@@H](C(N4[C@@H](C(NC(CSC1C)C(=O)NC(C(=O)NC(C(=O)NC(CC(NC(C(C2)C)=O)=O)C(=O)O)C)CSC3)=O)CCC4)=O)CC(C)C)=C)C(NC(=O)C=5N=C(C(NC(=O)C)C)SC5)CC=6NC=NC6',
                                      'name': 'Nocathioamide A',
                                      'reason': 'Detected 0 peptide bonds, '
                                                'need at least 5 for a '
                                                'polypeptide'}],
    'attempt': 2,
    'reasoning': 'To improve the classification of polypeptides, we need to '
                 'address a few key aspects from the previous implementation:\n'
                 '\n'
                 '1. **Peptide Bond Patterns**: While the existing pattern for '
                 'peptide bonds is a good start, it may not cover all '
                 'configurations of peptide bonds involving different side '
                 'chains. We should expand this to recognize not only linear '
                 'chains but also cycles and other configurations.\n'
                 '\n'
                 '2. **Backbone Pattern Recognition**: The recognition of the '
                 'characteristic polypeptide backbone should include more '
                 'flexible SMARTS patterns to accommodate diverse polypeptide '
                 'structures.\n'
                 '\n'
                 '3. **Molecular Complexity and Weight**: The criteria for '
                 'molecular weight and rotatable bonds may be too strict or '
                 'not sufficiently descriptive. Adjusting these thresholds or '
                 'using other descriptors (such as atom counts, specific atom '
                 'types like nitrogen for peptide groups) could improve '
                 'accuracy.\n'
                 '\n'
                 '4. **Handling False Negatives**: Review of the false '
                 'negatives shows that "Molecular complexity too low for a '
                 'polypeptide" was a common reason. This suggests that our '
                 'thresholds might be set too high or that our SMARTS patterns '
                 "aren't capturing the full variety of polypeptide "
                 'structures.\n'
                 '\n'
                 '5. **Full Connectivity Check**: Ensure the identified '
                 'peptide bonds create fully connected chains as expected for '
                 'polypeptides.\n'
                 '\n'
                 "Let's adjust the code taking these considerations into "
                 'account:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 251,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.7150997150997151,
    'f1': 0.8338870431893688,
    'accuracy': 0.7150997150997151,
    'negative_predictive_value': 0.0}