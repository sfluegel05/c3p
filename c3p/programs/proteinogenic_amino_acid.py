"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool, str: True and reason if molecule is a proteinogenic amino acid, otherwise False and reason
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for standard amino acid backbone (N-C-C(=O)-O) with chirality
    backbone_pattern = Chem.MolFromSmarts("N[C@H](C(=O)O)")

    # Check if the mol has the standard amino acid backbone with chirality
    if not mol.HasSubstructMatch(backbone_pattern):
        # Glycine exception, which is achiral
        glycine_pattern = Chem.MolFromSmarts("NCC(=O)O")
        if mol.HasSubstructMatch(glycine_pattern):
            return True, "SMILES matches glycine structure, an exception to chirality"
        return False, "No standard amino acid backbone with chirality found"

    # Check for chirality
    if not any(atom.GetChiralTag() != Chem.CHI_UNSPECIFIED for atom in mol.GetAtoms()):
        return False, "Chirality not found in molecule except for glycine"

    return True, "Valid backbone and chirality detected; matches a proteinogenic amino acid"

# Example test case
smiles_example = "N[C@@H](CC(N)=O)C(O)=O"  # L-asparagine
result, reason = is_proteinogenic_amino_acid(smiles_example)
print(f"Is proteinogenic amino acid: {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:83813',
                          'name': 'proteinogenic amino acid',
                          'definition': 'Any of the 23 alpha-amino acids that '
                                        'are precursors to proteins, and are '
                                        'incorporated into proteins during '
                                        'translation. The group includes the '
                                        '20 amino acids encoded by the nuclear '
                                        'genes of eukaryotes together with '
                                        'selenocysteine, pyrrolysine, and '
                                        'N-formylmethionine. Apart from '
                                        'glycine, which is non-chiral, all '
                                        'have L configuration.',
                          'parents': ['CHEBI:33709'],
                          'xrefs': ['Wikipedia:Proteinogenic_amino_acid'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Error: Python argument types in\n'
               '    Mol.HasSubstructMatch(Mol, NoneType)\n'
               'did not match C++ signature:\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle '
               'query, RDKit::SubstructMatchParameters params=True)\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
               'RDKit::SubstructMatchParameters params)\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle '
               'query, bool recursionPossible=True, bool useChirality=False, '
               'bool useQueryQueryMatches=False)\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
               'bool recursionPossible=True, bool useChirality=False, bool '
               'useQueryQueryMatches=False)\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'No standard amino acid '
                                               'backbone with chirality found'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'No standard amino acid '
                                               'backbone with chirality found'},
                                 {   'smiles': 'C(C(O)=O)C/C=C\\C/C=C\\C\\C=C/C=C/C=C/[C@H]1[C@H](C/C=C\\CC)O1',
                                     'name': '(16S,17S)-epoxy-(4Z,7Z,10Z,12E,14E,19Z)-docosahexaenoic '
                                             'acid',
                                     'reason': 'No standard amino acid '
                                               'backbone with chirality found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'No standard amino acid '
                                               'backbone with chirality found'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'No standard amino acid '
                                               'backbone with chirality found'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'No standard amino acid '
                                               'backbone with chirality found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'No standard amino acid '
                                               'backbone with chirality found'},
                                 {   'smiles': 'O=C1OC(O)C(=C1C(O)C(C)C)C',
                                     'name': '5-Hydroxy-3-(1-hydroxy-2-methylpropyl)-4-methyl-2(5H)-furanone',
                                     'reason': 'No standard amino acid '
                                               'backbone with chirality found'},
                                 {   'smiles': '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1',
                                     'name': '(+)-tephrorin B',
                                     'reason': 'No standard amino acid '
                                               'backbone with chirality found'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)CC4=CC=CC(=C4)C(=O)O)C5=CC=CC=C5N2C)[C@@H](C)CO',
                                     'name': 'LSM-9341',
                                     'reason': 'No standard amino acid '
                                               'backbone with chirality '
                                               'found'}],
    'sample_false_negatives': [],
    'attempt': 2,
    'reasoning': 'To address the problem and improve the classification of '
                 "proteinogenic amino acids using their SMILES strings, here's "
                 'a step-by-step analysis and potential solution:\n'
                 '\n'
                 '1. **Understand the Problem**: The goal is to identify '
                 'whether a given molecule is a proteinogenic amino acid. This '
                 'includes checking for the presence of a standard peptide '
                 'backbone structure, recognizing chirality for most amino '
                 'acids, and accounting for the special cases such as '
                 'glycine.\n'
                 '\n'
                 '2. **Identify Issues in Previous Code**:\n'
                 '   - The error indicates that a `NoneType` was used in a '
                 'substructure match, suggesting that the SMARTS pattern '
                 "didn't compile correctly or wasn't used correctly.\n"
                 '   - Ensure that chirality and substructure checks are '
                 'correctly implemented.\n'
                 '\n'
                 '3. **Chirality Considerations**: Most amino acids are chiral '
                 'with an L-configuration; however, glycine is an exception as '
                 "it is non-chiral. We'll need to handle this by checking "
                 'chirality carefully.\n'
                 '\n'
                 '4. **Revised Strategy**:\n'
                 '   - Validate the SMILES string using RDKit.\n'
                 '   - Define and correctly use the SMARTS patterns for both '
                 'the general peptide backbone and special cases.\n'
                 '   - Check for chirality using RDKit features, ensuring '
                 'non-chiral exceptions are handled.\n'
                 '   - Provide detailed reasons for the classification results '
                 'for better understanding.\n'
                 '\n'
                 'Here is a revised Python program that addresses these '
                 'considerations:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 27,
    'num_false_positives': 100,
    'num_true_negatives': 1017,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.2125984251968504,
    'recall': 1.0,
    'f1': 0.35064935064935066,
    'accuracy': 0.9125874125874126,
    'negative_predictive_value': 1.0}