"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
"""
Classifies: CHEBI:85259 alpha-amino acid ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester is formed by the condensation of an alpha-amino acid with an alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the ester group (-COO-)
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) == 0:
        return False, "No ester group found"

    # Look for the alpha-amino acid pattern (NH2-CH-COOH)
    amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2][CX4H][CX3](=[OX1])[OX2H1]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if len(amino_acid_matches) == 0:
        return False, "No alpha-amino acid moiety found"

    # Check if the ester group is connected to the alpha-carbon of the amino acid
    alpha_carbon = None
    for match in amino_acid_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'C' and atom.GetDegree() == 4:  # Alpha carbon
                alpha_carbon = atom_idx
                break
        if alpha_carbon is not None:
            break

    if alpha_carbon is None:
        return False, "No alpha-carbon found in the amino acid moiety"

    # Check if the ester group is connected to the alpha-carbon
    ester_connected = False
    for match in ester_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'C' and atom.GetIdx() == alpha_carbon:
                ester_connected = True
                break
        if ester_connected:
            break

    if not ester_connected:
        return False, "Ester group not connected to the alpha-carbon of the amino acid"

    return True, "Contains alpha-amino acid moiety with ester group connected to the alpha-carbon"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:85259',
                          'name': 'alpha-amino acid ester',
                          'definition': 'The amino acid ester derivative obtained the formal condensation of an alpha-amino acid with an alcohol.',
                          'parents': ['CHEBI:85258', 'CHEBI:85260']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}