"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: CHEBI:51294 N-acylsphingosine
"""
from rdkit import Chem

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    An N-acylsphingosine is a sphingosine backbone with an acyl group attached to the nitrogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphingosine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for N-acylsphingosine
    n_acyl_sphingosine_smarts = "[NX3;!$(N-[N])]-C(=O)-[C@H]([O])[C@H]([O])[*]"
    n_acyl_sphingosine_pattern = Chem.MolFromSmarts(n_acyl_sphingosine_smarts)
    if n_acyl_sphingosine_pattern is None:
        return False, "Invalid SMARTS pattern"

    if not mol.HasSubstructMatch(n_acyl_sphingosine_pattern):
        return False, "Molecule does not match N-acylsphingosine pattern"

    # Check for long hydrocarbon chain (minimum 12 carbons)
    # Starting from the beta carbon
    matches = mol.GetSubstructMatches(n_acyl_sphingosine_pattern)
    for match in matches:
        # Index positions in the SMARTS pattern
        # 0 - Nitrogen
        # 1 - Carbonyl carbon
        # 2 - Alpha carbon [C@H]([O])
        # 3 - Beta carbon [C@H]([O])
        # 4 - Next atom (could be carbon)
        beta_carbon_idx = match[3]
        beta_carbon = mol.GetAtomWithIdx(beta_carbon_idx)

        # Traverse the chain from beta carbon
        visited = set()
        chain_length = 0
        atoms_to_visit = [beta_carbon]
        while atoms_to_visit:
            atom = atoms_to_visit.pop()
            idx = atom.GetIdx()
            if idx in visited or idx in match:
                continue
            visited.add(idx)
            if atom.GetAtomicNum() == 6:
                chain_length += 1
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:
                        atoms_to_visit.append(neighbor)
        if chain_length >= 12:
            return True, "Molecule is an N-acylsphingosine with sphingosine backbone"
        else:
            return False, f"Hydrocarbon chain too short (length {chain_length})"

    return False, "Does not match N-acylsphingosine structure"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:51294',
        'name': 'N-acylsphingosine',
        'definition': 'The parent compounds of the ceramide family, composed of '
                      'sphingosine having an unspecified fatty acyl group attached '
                      'to the nitrogen.',
        'parents': ['CHEBI:76165', 'CHEBI:76158']
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
    'accuracy': None
}