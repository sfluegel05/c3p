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

    # Look for the alpha-amino acid pattern (NH2/NH/N-CH-COOH/COOR)
    # More flexible pattern that accounts for substituted amines and different carboxylate forms
    amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0][CX4H][CX3](=[OX1])[OX2]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if len(amino_acid_matches) == 0:
        return False, "No alpha-amino acid moiety found"

    # Check if the ester group is connected to the carboxylate of the amino acid
    ester_connected = False
    for aa_match in amino_acid_matches:
        # Get the carboxylate carbon and oxygen indices
        carboxyl_c = aa_match[2]
        carboxyl_o = aa_match[3]
        
        # Check if any ester group is connected to the carboxylate oxygen
        for ester_match in ester_matches:
            ester_o = ester_match[1]
            if ester_o == carboxyl_o:
                ester_connected = True
                break
        if ester_connected:
            break

    if not ester_connected:
        return False, "Ester group not connected to the carboxylate of the amino acid"

    # Check if the amino acid has an alpha-carbon with at least one hydrogen
    alpha_carbon_valid = False
    for aa_match in amino_acid_matches:
        alpha_c = aa_match[1]
        atom = mol.GetAtomWithIdx(alpha_c)
        if atom.GetTotalNumHs() > 0:
            alpha_carbon_valid = True
            break

    if not alpha_carbon_valid:
        return False, "No valid alpha-carbon found in the amino acid moiety"

    return True, "Contains alpha-amino acid moiety with ester group connected to the carboxylate"


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