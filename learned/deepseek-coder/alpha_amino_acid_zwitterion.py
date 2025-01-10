"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
"""
Classifies: CHEBI:57925 alpha-amino-acid zwitterion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    An alpha-amino-acid zwitterion has a protonated amino group ([NH3+]) and a deprotonated carboxyl group ([O-])
    attached to the same alpha-carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino-acid zwitterion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the zwitterion pattern: [NH3+] attached to a carbon, which is also attached to [O-]
    zwitterion_pattern = Chem.MolFromSmarts("[NH3+][CX4][CX3](=[OX1])[OX1-]")
    if not mol.HasSubstructMatch(zwitterion_pattern):
        return False, "No zwitterion pattern found ([NH3+]-C-C(=O)[O-])"

    # Verify that the central carbon (alpha-carbon) is attached to both [NH3+] and [O-]
    alpha_carbon_pattern = Chem.MolFromSmarts("[CX4]([NH3+])([CX3](=[OX1])[OX1-])")
    if not mol.HasSubstructMatch(alpha_carbon_pattern):
        return False, "No alpha-carbon with both [NH3+] and [O-] groups found"

    # Check that the alpha-carbon has at least 3 bonds (including H)
    alpha_carbon_matches = mol.GetSubstructMatches(alpha_carbon_pattern)
    for match in alpha_carbon_matches:
        alpha_carbon_idx = match[0]  # Index of the alpha-carbon
        alpha_carbon_atom = mol.GetAtomWithIdx(alpha_carbon_idx)
        if alpha_carbon_atom.GetDegree() < 3:  # Alpha-carbon should have at least 3 bonds
            return False, "Alpha-carbon has too few bonds"

    # Remove molecular weight restriction as it's not a defining characteristic

    return True, "Contains [NH3+]-C-C(=O)[O-] pattern with alpha-carbon in zwitterionic form"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:57925',
        'name': 'alpha-amino-acid zwitterion',
        'definition': 'An amino acid-zwitterion obtained by transfer of a proton from the carboxy to the amino group of any alpha-amino acid; major species at pH 7.3.',
        'parents': ['CHEBI:33704', 'CHEBI:57924']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
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
    'accuracy': 0.9998521228585199
}