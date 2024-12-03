"""
Classifies: CHEBI:60980 beta-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-glucoside.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-glucoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for beta-D-glucopyranose structure
    beta_glucose_smiles = 'O[C@@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O1)CO'
    beta_glucose_mol = Chem.MolFromSmiles(beta_glucose_smiles)

    if beta_glucose_mol is None:
        return False, "Failed to create beta-glucose molecule"

    # Perform substructure search
    if not mol.HasSubstructMatch(beta_glucose_mol):
        return False, "No beta-D-glucopyranose structure found"

    # Check the anomeric carbon for beta configuration
    matches = mol.GetSubstructMatches(beta_glucose_mol)
    for match in matches:
        anomeric_carbon_idx = match[0]  # The first atom in the match is the anomeric carbon
        anomeric_carbon = mol.GetAtomWithIdx(anomeric_carbon_idx)
        if anomeric_carbon.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
            return True, "Beta configuration confirmed at the anomeric carbon"

    return False, "Beta configuration not confirmed at the anomeric carbon"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:60980',
                          'name': 'beta-glucoside',
                          'definition': 'A glucoside in which the anomeric '
                                        'carbon of the glycosidic bond is in a '
                                        'beta configuration',
                          'parents': ['CHEBI:24278']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 75,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}