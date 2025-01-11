"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: trans-2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    These are CoA thioesters with a trans double bond at position 2 of the acyl chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trans-2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Essential structural patterns
        patterns = {
            # Core CoA structure components
            "adenine": Chem.MolFromSmarts("c1nc(N)c2ncnc2n1"),
            "phosphate": Chem.MolFromSmarts("OP(=O)(O)O"),
            "thioester": Chem.MolFromSmarts("C(=O)SC"),
            "pantetheine": Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)"),
            # Pattern for alpha-beta unsaturated thioester, without stereochemistry
            "enoyl": Chem.MolFromSmarts("C=CC(=O)SC")
        }

        # Check for CoA moiety components
        if not mol.HasSubstructMatch(patterns["adenine"]):
            return False, "No CoA moiety found (missing adenine)"

        phosphate_matches = len(mol.GetSubstructMatches(patterns["phosphate"]))
        if phosphate_matches < 2:
            return False, "No CoA moiety found (insufficient phosphate groups)"

        if not mol.HasSubstructMatch(patterns["thioester"]):
            return False, "No thioester linkage found"

        if not mol.HasSubstructMatch(patterns["pantetheine"]):
            return False, "Missing pantetheine arm of CoA"

        # Find all double bonds connected to thioester
        enoyl_matches = mol.GetSubstructMatches(patterns["enoyl"])
        if not enoyl_matches:
            return False, "No α,β-unsaturated thioester found"

        # Check stereochemistry of the double bonds
        found_trans = False
        for match in enoyl_matches:
            # Get the double bond atoms
            double_bond_atoms = [match[0], match[1]]  # C=C atoms from the match
            bond = mol.GetBondBetweenAtoms(double_bond_atoms[0], double_bond_atoms[1])
            
            # Check if it's a double bond and has stereochemistry
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                # E/trans configuration is represented by STEREOCIS in RDKit
                if bond.GetStereo() == Chem.BondStereo.STEREOE:
                    found_trans = True
                    break

        if not found_trans:
            return False, "No trans configuration found at position 2"

        # Additional check for position (should be alpha to the thioester)
        # This is already enforced by our enoyl pattern

        return True, "Contains CoA moiety with trans-2-enoyl group"

    except Exception as e:
        return False, f"Error in structure analysis: {str(e)}"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50998',
                          'name': 'trans-2-enoyl-CoA',
                          'definition': 'An unsaturated fatty acyl-CoA that '
                                        'results from the formal condensation '
                                        'of the thiol group of coenzyme A with '
                                        'the carboxy group of any '
                                        '2,3-trans-enoic acid.',
                          'parents': ['CHEBI:19573'],
                          'xrefs': ['KEGG:C00658'],
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
               'CCCCCCCC\\C=C/CCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (2E,9Z)-octadecadienoyl-CoA REASON: MISSED Error in '
               'structure analysis: Invalid SMARTS pattern: \\C=C\\C(=O)S\n'
               ' * SMILES: '
               'CCCCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: trans-dodec-2-enoyl-CoA REASON: MISSED Error in '
               'structure analysis: Invalid SMARTS pattern: \\C=C\\C(=O)S\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: trans-2-triacontenoyl-CoA REASON: MISSED Error in '
               'structure analysis: Invalid SMARTS pattern: \\C=C\\C(=O)S\n'
               ' * SMILES: '
               'CCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: trans-dec-2-enoyl-CoA REASON: MISSED Error in structure '
               'analysis: Invalid SMARTS pattern: \\C=C\\C(=O)S\n'
               ' * SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)/C=C/CCCCCCCCCCCCC)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O '
               'NAME: (E)-hexadec-2-enoyl-CoA REASON: MISSED Error in '
               'structure analysis: Invalid SMARTS pattern: \\C=C\\C(=O)S\n'
               ' * SMILES: '
               'CCCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: trans-2-undecenoyl-CoA REASON: MISSED Error in structure '
               'analysis: Invalid SMARTS pattern: \\C=C\\C(=O)S\n'
               ' * SMILES: '
               'CCCCCC\\C=C/CCCCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (2E,13Z)-icosadienoyl-CoA REASON: MISSED Error in '
               'structure analysis: Invalid SMARTS pattern: \\C=C\\C(=O)S\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: trans-2-tetracosenoyl-CoA REASON: MISSED Error in '
               'structure analysis: Invalid SMARTS pattern: \\C=C\\C(=O)S\n'
               ' * SMILES: '
               'CCCCCCCC\\C=C\\CCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (2E,9E)-octadecadienoyl-CoA REASON: MISSED Error in '
               'structure analysis: Invalid SMARTS pattern: \\C=C\\C(=O)S\n'
               ' * SMILES: '
               'CCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: trans-hex-2-enoyl-CoA REASON: MISSED Error in structure '
               'analysis: Invalid SMARTS pattern: \\C=C\\C(=O)S\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C1C=C(N=C2N1NN=N2)C3=CC=CC=C3',
                                     'name': '5-phenyl-1,7-dihydrotetrazolo[1,5-a]pyrimidine',
                                     'reason': 'No CoA moiety found (missing '
                                               'adenine)'},
                                 {   'smiles': 'O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)[C@H](COC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\CCCC',
                                     'name': 'TG(14:1(9Z)/22:1(13Z)/18:1(11Z))',
                                     'reason': 'No CoA moiety found (missing '
                                               'adenine)'},
                                 {   'smiles': 'CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H]1O',
                                     'name': 'beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-{alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->6)]-alpha-D-Manp-(1->6)}-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)-beta-D-GlcpNAc',
                                     'reason': 'No CoA moiety found (missing '
                                               'adenine)'},
                                 {   'smiles': '[C@@]12([C@@]([C@]([C@@H](CC1)C)(CCC(CC)C)C)(CCC[C@@H]2C)[H])C',
                                     'name': 'clerodane',
                                     'reason': 'No CoA moiety found (missing '
                                               'adenine)'},
                                 {   'smiles': 'S(O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]([C@@H](O)[C@H](O)CO[C@]2(O[C@H]([C@H](NC(=O)C)[C@@H](O)C2)[C@H](O)[C@H](O)CO)C(O)=O)[C@@H](NC(=O)C)CO)O[C@@H]1CO)(O)(=O)=O',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4R,5S)-5-acetamido-4-[(2R,3R,4R,5R,6R)-3,4-dihydroxy-6-(hydroxymethyl)-5-sulfooxyoxan-2-yl]oxy-2,3,6-trihydroxyhexoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No CoA moiety found (missing '
                                               'adenine)'},
                                 {   'smiles': 'O=C(NC(CC(C)C)C(=O)N[C@@H](CCCN=C(N)N)C=O)C(NC(=O)CC)CC(C)C',
                                     'name': 'Leupeptin Pr-LL',
                                     'reason': 'No CoA moiety found (missing '
                                               'adenine)'},
                                 {   'smiles': 'COC1=CC=C(C=C1)CCNC(=O)C(=C2C3=CC=CC=C3C(=N2)NC(=O)C4=CC=CS4)C#N',
                                     'name': 'N-[3-[1-cyano-2-[2-(4-methoxyphenyl)ethylamino]-2-oxoethylidene]-1-isoindolyl]-2-thiophenecarboxamide',
                                     'reason': 'No CoA moiety found (missing '
                                               'adenine)'},
                                 {   'smiles': 'O=C1O[C@@H](CC=2C1=C(OC)C(O)=CC2)CCC[C@@H](O)C',
                                     'name': 'Penicimarin C',
                                     'reason': 'No CoA moiety found (missing '
                                               'adenine)'},
                                 {   'smiles': 'CO[C@@H]1[C@@H]2[C@@H](C[C@@H]3[C@@](O2)([C@H]([C@@H]([C@H](O3)C(=O)OC)O)O)O)OC1N4C=NC5=C4N=CN=C5N',
                                     'name': 'LSM-4497',
                                     'reason': 'No CoA moiety found '
                                               '(insufficient phosphate '
                                               'groups)'},
                                 {   'smiles': 'C[C@@H]1[C@H](O)CCC2=CC[C@H](C[C@]12C)C(C)=C',
                                     'name': '1-deoxycapsidiol',
                                     'reason': 'No CoA moiety found (missing '
                                               'adenine)'}],
    'sample_false_negatives': [],
    'attempt': 2,
    'reasoning': "Let's analyze what went wrong and improve the code:\n"
                 '\n'
                 '1. Main Issue: The error messages show that the SMARTS '
                 'patterns for trans/cis double bonds were invalid. The '
                 'backslashes in the SMARTS patterns were causing problems '
                 'because they need to be properly escaped in Python strings.\n'
                 '\n'
                 '2. Improvements needed:\n'
                 '- Fix the SMARTS patterns for trans/cis configurations by '
                 'using raw strings (r"...")\n'
                 '- Simplify the double bond detection approach\n'
                 '- Make the pattern matching more robust\n'
                 '- Add better handling of stereochemistry\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 26,
    'num_false_positives': 100,
    'num_true_negatives': 67769,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.20634920634920634,
    'recall': 1.0,
    'f1': 0.34210526315789475,
    'accuracy': 0.9985271374917152,
    'negative_predictive_value': 1.0}