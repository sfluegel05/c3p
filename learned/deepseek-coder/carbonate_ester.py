"""
Classifies: CHEBI:46722 carbonate ester
"""
"""
Classifies: CHEBI:51052 carbonate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester is carbonic acid in which the hydrogens have been replaced by organyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the carbonate ester functional group pattern: [O]-[C](=O)-[O]
    carbonate_ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[OX2]")
    if not mol.HasSubstructMatch(carbonate_ester_pattern):
        return False, "No carbonate ester functional group found"

    # Check that the oxygen atoms in the carbonate ester group are bonded to carbon atoms (organyl groups)
    for match in mol.GetSubstructMatches(carbonate_ester_pattern):
        # Get the atoms in the match
        atom1 = mol.GetAtomWithIdx(match[0])  # First oxygen
        atom2 = mol.GetAtomWithIdx(match[2])  # Second oxygen

        # Check that both oxygens are bonded to carbon atoms (organyl groups)
        # Iterate through the neighbors of each oxygen atom
        for oxygen_atom in [atom1, atom2]:
            has_carbon_neighbor = False
            for neighbor in oxygen_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon atom
                    has_carbon_neighbor = True
                    break
            if not has_carbon_neighbor:
                return False, f"Oxygen atom {oxygen_atom.GetIdx()} is not bonded to a carbon atom (not an organyl group)"

    # Additional check to avoid false positives: ensure no other interfering functional groups
    # For example, check for carboxylic acids, esters, etc.
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Contains a carboxylic acid group, which is not a carbonate ester"

    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) > 1:
        return False, "Contains multiple ester groups, which may interfere with carbonate ester classification"

    return True, "Contains a carbonate ester functional group with organyl substituents"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:51052',
        'name': 'carbonate ester',
        'definition': 'Any carbonate that is carbonic acid in which the hydrogens have been replaced by organyl groups.',
        'parents': ['CHEBI:51050', 'CHEBI:51051']
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.0
}