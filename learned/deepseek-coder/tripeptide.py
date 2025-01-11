"""
Classifies: CHEBI:47923 tripeptide
"""
"""
Classifies: CHEBI:48030 tripeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide consists of three amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tripeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for peptide bonds (amide linkages)
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3][CX4]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) != 2:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, need exactly 2"

    # Check for N-terminal (NH2 or H)
    n_terminal_pattern = Chem.MolFromSmarts("[NX3;H2,H1]")
    if not mol.HasSubstructMatch(n_terminal_pattern):
        return False, "No N-terminal group found"

    # Check for C-terminal (COOH or O)
    c_terminal_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(c_terminal_pattern):
        return False, "No C-terminal group found"

    # Count amino acid residues (alpha carbons with NH and CO groups)
    alpha_carbon_pattern = Chem.MolFromSmarts("[CX4H]([NX3H1,CX3](=O))[CX4H1]")
    alpha_carbon_matches = mol.GetSubstructMatches(alpha_carbon_pattern)
    if len(alpha_carbon_matches) != 3:
        return False, f"Found {len(alpha_carbon_matches)} amino acid residues, need exactly 3"

    return True, "Contains three amino acid residues connected by two peptide bonds"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:48030',
                          'name': 'tripeptide',
                          'definition': 'Any oligopeptide that consists of three '
                                        'amino-acid residues connected by peptide '
                                        'linkages.',
                          'parents': ['CHEBI:25676', 'CHEBI:26561']},
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