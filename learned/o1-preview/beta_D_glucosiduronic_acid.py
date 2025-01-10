"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
"""
Classifies: beta-D-glucosiduronic acid derivatives
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid derivative based on its SMILES string.
    A beta-D-glucosiduronic acid derivative contains a beta-D-glucuronic acid moiety attached via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronic acid derivative, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for beta-D-glucuronic acid moiety attached via glycosidic bond
    pattern = Chem.MolFromSmarts("""
    [C@@H]1([O][#6])O[C@@H]([C@@H](O)[C@H](O)[C@H](O)[C@H]1O)C(=O)[O]
    """)
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check for beta-D-glucuronic acid moiety
    if not mol.HasSubstructMatch(pattern):
        return False, "Beta-D-glucuronic acid moiety not found"

    # Check if the anomeric carbon (C1) is involved in a glycosidic bond
    # Find matches to the pattern
    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "Beta-D-glucuronic acid moiety not found"

    # Check each match to see if the glycosidic bond exists
    for match in matches:
        anomeric_carbon_idx = match[0]  # Index of the anomeric carbon (C1)
        anomeric_carbon = mol.GetAtomWithIdx(anomeric_carbon_idx)
        # Get bonds connected to the anomeric carbon
        bonds = anomeric_carbon.GetBonds()
        glycosidic_bond_found = False
        for bond in bonds:
            # Check if bond is to an oxygen atom that connects to another group
            neighbor = bond.GetOtherAtom(anomeric_carbon)
            if neighbor.GetAtomicNum() == 8:
                # Check if oxygen atom is connected to another heavy atom outside the glucuronic acid unit
                oxygen_atom = neighbor
                for obond in oxygen_atom.GetBonds():
                    other = obond.GetOtherAtom(oxygen_atom)
                    if other.GetIdx() != anomeric_carbon_idx and other.GetAtomicNum() > 1:
                        glycosidic_bond_found = True
                        break
                if glycosidic_bond_found:
                    break
        if glycosidic_bond_found:
            return True, "Beta-D-glucuronic acid moiety with glycosidic bond found"

    return False, "Beta-D-glucuronic acid moiety without glycosidic bond found"

__metadata__ = {   
    'chemical_class': {   
        'id': None,
        'name': 'beta-D-glucosiduronic acid',
        'definition': 'A glucosiduronic acid resulting from the formal condensation of any substance with beta-D-glucuronic acid to form a glycosidic bond.',
        'parents': []
        },
    'config': {   
        'llm_model_name': None,
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