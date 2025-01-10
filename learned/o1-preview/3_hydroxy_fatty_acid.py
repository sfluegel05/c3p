"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: 3-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid is a fatty acid with a hydroxy group on the third carbon from the carboxyl end.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure molecule is sanitized
    try:
        Chem.SanitizeMol(mol)
    except Chem.rdchem.KekulizeException:
        return False, "Molecule could not be sanitized"

    # Add explicit hydrogens
    mol = Chem.AddHs(mol)

    # Check for carboxylic acid group (allowing for protonated and deprotonated forms)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H0,H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for hydrocarbon chain (alkyl chain)
    alkyl_chain_pattern = Chem.MolFromSmarts("[CX4H2][CX4H2]")
    if not mol.HasSubstructMatch(alkyl_chain_pattern):
        return False, "No alkyl chain found"

    # Check for hydroxy group on the 3-position from the carboxyl group
    # Define a pattern where a carboxylic acid carbon is connected to the chain with a hydroxy on the 3rd carbon
    three_hydroxy_fatty_acid_pattern = Chem.MolFromSmarts("""
    [CX3](=O)[OX1H0,H1][CH2][CH](O)[CH2]
    """)

    # The above pattern matches a carboxylic acid connected to a chain where the third carbon has a hydroxy group

    if not mol.HasSubstructMatch(three_hydroxy_fatty_acid_pattern):
        return False, "No 3-hydroxy group found in fatty acid chain"

    # Optionally, check for long chain length typical of fatty acids
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 4:
        return False, f"Chain too short for fatty acid (found {num_carbons} carbons)"

    return True, "Molecule is a 3-hydroxy fatty acid with carboxylic acid group and hydroxy on 3-position"

__metadata__ = {   'chemical_class': {   'id': None,
                              'name': '3-hydroxy fatty acid',
                              'definition': 'Any fatty acid with a hydroxy functional group in the beta- or 3-position. Beta-hydroxy fatty acids accumulate during cardiac hypoxia and can also be used as chemical markers of bacterial endotoxins.',
                              'parents': ['fatty acid']},
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
        'attempt': 1,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None,
        'num_true_positives': None,
        'num_false_positives': None,
        'num_true_negatives': None,
        'num_false_negatives': None,
        'num_negatives': None,
        'precision': None,
        'recall': None,
        'f1': None,
        'accuracy': None}