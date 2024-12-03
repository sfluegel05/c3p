"""
Classifies: CHEBI:136889 5beta steroid
"""
from rdkit import Chem

def is_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 5beta steroid (steroid with beta-configuration at position 5).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 5beta steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the steroid core structure
    steroid_core = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3CCC(C4)C")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "Molecule does not have a steroid core structure"

    # Check for beta-configuration at position 5
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    position_5_chiral = None
    for center in chiral_centers:
        if center[0] == 4:  # Position 5 in the SMILES string corresponds to index 4 in RDKit
            position_5_chiral = center
            break

    if position_5_chiral is None:
        return False, "No chiral center found at position 5"

    if position_5_chiral[1] != 'S':
        return False, "Chiral center at position 5 is not in beta-configuration"

    return True, "Molecule is a 5beta steroid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:136889',
                          'name': '5beta steroid',
                          'definition': 'Any steroid that has '
                                        'beta-configuration at position 5.',
                          'parents': ['CHEBI:35341']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[19:19:57] SMILES Parse Error: syntax error while parsing: '
             'O1[C@@H](C[C@H]([C@@]2([C@@]3([C@@](CC2)(/C(/CCC3)=C/C=C\x04/C[C@@H](O)C[C@H](O)C4=C)[H])C)[H])C)C[C@@](O)(C1=O)C\n'
             '[19:19:57] SMILES Parse Error: Failed parsing SMILES '
             "'O1[C@@H](C[C@H]([C@@]2([C@@]3([C@@](CC2)(/C(/CCC3)=C/C=C\x04/C[C@@H](O)C[C@H](O)C4=C)[H])C)[H])C)C[C@@](O)(C1=O)C' "
             'for input: '
             "'O1[C@@H](C[C@H]([C@@]2([C@@]3([C@@](CC2)(/C(/CCC3)=C/C=C\x04/C[C@@H](O)C[C@H](O)C4=C)[H])C)[H])C)C[C@@](O)(C1=O)C'\n",
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 48,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}