"""
Classifies: CHEBI:25676 oligopeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is a peptide containing a relatively small number of amino acids (usually <20).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligopeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for amide bonds (-C(=O)NH- patterns)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    # Count number of amide bonds to estimate chain length
    num_amides = len(amide_matches)
    
    # Each amide bond typically corresponds to one amino acid in the chain (except for the terminal residues)
    num_amino_acids = num_amides + 1  # Assuming a linear chain, last amino acid may not form an additional amide bond
    
    # Threshold for oligopeptides: typically less than 20 amino acids
    if num_amino_acids < 20:
        return True, f"Contains {num_amino_acids} amino acids, which is typical of an oligopeptide"
    else:
        return False, f"Contains {num_amino_acids} amino acids, exceeding typical oligopeptide range"

    return False, "Undefined classification logic"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25676',
                          'name': 'oligopeptide',
                          'definition': 'A peptide containing a relatively '
                                        'small number of amino acids.',
                          'parents': ['CHEBI:16670'],
                          'xrefs': ['Wikipedia:Oligopeptide'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O=C(N1[C@H](C(=O)N[C@H](C(=O)NC(C(=O)NC(C(=O)N[C@H](C(=O)NC(C(=O)NC(CNCCO)CC2=CC=CC=C2)CCC(=O)N)CCC(=O)O)(C)C)(C)C)C(C)C)CCC1)C(NC(=O)[C@@H](NC(=O)CNC(=O)C(NC(=O)[C@@H](NC(=O)C(NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)C(NC(=O)[C@@H](NC(=O)C(NC(=O)[C@@H]3N(C(=O)C(NC(=O)C)(C)C)CCC3)(C)C)C)(C)C)C)CCC(=O)N)(C)C)C(C)C)(C)C)CC(C)C)(C)C',
                                     'name': 'Septocylindrin A',
                                     'reason': 'Contains 23 amino acids, '
                                               'exceeding typical oligopeptide '
                                               'range'}],
    'sample_false_negatives': [   {   'smiles': 'C(N[C@H](C(N[C@@H](CC(=O)N)C(N1CCC[C@H]1C(=O)N[C@H](C(=O)N[C@H](C(=O)N2[C@H](C(=O)N[C@H](C(N[C@H](C(N[C@@H](CC(=O)N)C(=O)N3[C@H](C(=O)N[C@H](C(N[C@H](C(N4[C@H](C(=O)N[C@H](C(N[C@H](C(N[C@H](C(=O)N)CCCCN)=O)CCC(=O)N)=O)CO)CCC4)=O)CO)=O)CC=5C=CC(=CC5)O)CCC3)=O)CCC(=O)N)=O)CO)CCC2)CO)CC=6C=CC(=CC6)O)=O)=O)CCCNC(=N)N)([C@H](CCC(=O)O)NC(=O)C)=O',
                                      'name': 'Ac-Glu-Arg-Asp-Pro-Tyr-Ser-Pro-Ser-Gln-Asp-Pro-Tyr-Ser-Pro-Ser-Gln-Lys-NH2',
                                      'reason': 'Contains 23 amino acids, '
                                                'exceeding typical '
                                                'oligopeptide range'},
                                  {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)N2C(C(=O)NC(C(=O)NC(C(=O)O)CO)CC3=CC=CC=C3)CCC2)(C)C)CC(C1)O)C(NC(=O)C(NC(=O)C4N(C(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC5=CC=CC=C5)(CC)C)CCC(=O)N)(C)C)C(CC)C)C(O)C)(C)C)CC(C)C)(C)C)CCC4)CCC(=O)N)(C)C',
                                      'name': 'Cephaibol Q',
                                      'reason': 'Contains 20 amino acids, '
                                                'exceeding typical '
                                                'oligopeptide range'},
                                  {   'smiles': 'C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(=O)N[C@H](C(N[C@H](C(N[C@H](C(N[C@H](C(N[C@@H](CC(=O)N)C(=O)N[C@H](C(N1CCC[C@H]1C(=O)N[C@H](C(N[C@H](C(N[C@H](C(=O)N)CCCCN)=O)CCC(=O)N)=O)CCC(=O)N)=O)CC(C)C)=O)CCCNC(=N)N)=O)CC(C)C)=O)CCC(=O)O)=O)CCCNC(=N)N)CCCCN)=O)CC=2C=CC=CC2)=O)CCC(=O)N)=O)CCC(=O)N)=O)CCC(=O)O)=O)CCC(=O)N)([C@H](CCC(=O)N)NC(=O)C)=O',
                                      'name': 'Ac-Gln-Gln-Glu-Gln-Gln-Phe-Lys-Arg-Glu-Leu-Arg-Asn-Leu-Pro-Gln-Gln-Lys-NH2',
                                      'reason': 'Contains 26 amino acids, '
                                                'exceeding typical '
                                                'oligopeptide range'},
                                  {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)N2C(C(=O)NC(C(=O)NC(C(=O)O)CO)CC3=CC=CC=C3)CCC2)(C)C)CC(C1)O)C(NC(=O)C(NC(=O)C4N(C(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC5=CC=CC=C5)(CC)C)CCC(=O)N)(C)C)C(CC)C)C(O)C)(C)C)CC(C)C)(C)C)CC(C4)O)CCC(=O)N)(C)C',
                                      'name': 'Cephaibol P',
                                      'reason': 'Contains 20 amino acids, '
                                                'exceeding typical '
                                                'oligopeptide range'}],
    'attempt': 0,
    'reasoning': 'To classify oligopeptides based on their SMILES strings, we '
                 'need to identify key features that are characteristic of '
                 'peptides and then determine if they meet the criteria of '
                 'being an "oligopeptide," which is typically defined as a '
                 'peptide consisting of a relatively small number of amino '
                 'acids (commonly less than 20).\n'
                 '\n'
                 '### Key Features of Oligopeptides:\n'
                 '1. **Amide Bonds**: Oligopeptides are sequences of amino '
                 'acids linked by amide (peptide) bonds.\n'
                 '2. **Amino Acids**: These are specific building blocks, each '
                 'containing an amine group and a carboxylic acid group.\n'
                 '3. **Chain Length**: An oligopeptide is generally shorter, '
                 'often containing less than 20-25 amino acids.\n'
                 '\n'
                 "We'll use RDKit to analyze SMILES strings and identify these "
                 "characteristics. We'll check for a series of amide bonds and "
                 'ensure the chain length is within the oligopeptide range. \n'
                 '\n'
                 "For this task, we'll create a function `is_oligopeptide`.",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 3591,
    'num_false_positives': 100,
    'num_true_negatives': 1,
    'num_false_negatives': 4,
    'num_negatives': None,
    'precision': 0.9729070712544026,
    'recall': 0.9988873435326843,
    'f1': 0.9857260499588251,
    'accuracy': 0.9718614718614719,
    'negative_predictive_value': 0.2}