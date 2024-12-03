"""
Classifies: CHEBI:23899 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains 20 carbon atoms in a chain
    carbon_chain_lengths = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            carbon_chain_length = 0
            visited = set()
            stack = [(atom, 0)]
            while stack:
                current_atom, length = stack.pop()
                if current_atom.GetIdx() in visited:
                    continue
                visited.add(current_atom.GetIdx())
                if current_atom.GetSymbol() == 'C':
                    carbon_chain_length = max(carbon_chain_length, length + 1)
                for neighbor in current_atom.GetNeighbors():
                    stack.append((neighbor, length + 1))
            carbon_chain_lengths.append(carbon_chain_length)

    if not any(length >= 20 for length in carbon_chain_lengths):
        return False, "No 20-carbon chain found"

    # Check for the presence of functional groups indicative of oxidation (e.g., hydroxyl, carbonyl, carboxyl)
    functional_groups = {
        'hydroxyl': Chem.MolFromSmarts('[OX2H]'),
        'carbonyl': Chem.MolFromSmarts('[CX3](=O)'),
        'carboxyl': Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    }

    for group_name, group_smarts in functional_groups.items():
        if mol.HasSubstructMatch(group_smarts):
            return True, f"Contains 20-carbon chain and {group_name} group"

    return False, "Does not contain necessary functional groups for icosanoids"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23899',
                          'name': 'icosanoid',
                          'definition': 'Any member of the group of signalling '
                                        'molecules arising from oxidation of '
                                        'the three C20 essential fatty acids '
                                        '(EFAs) icosapentaenoic acid (EPA), '
                                        'arachidonic acid (AA) and '
                                        'dihomo-gamma-linolenic acid (DGLA).',
                          'parents': ['CHEBI:61697']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 67,
    'num_false_positives': 12,
    'num_true_negatives': 8,
    'num_false_negatives': 25,
    'precision': 0.8481012658227848,
    'recall': 0.7282608695652174,
    'f1': 0.7836257309941521,
    'accuracy': None}