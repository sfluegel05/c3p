"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: CHEBI:33504 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide is a nucleoside phosphate resulting from the condensation of the 3' or 5' hydroxy group of a nucleoside with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for nucleobases
    purine_pattern = Chem.MolFromSmarts('c1ncnc2nccc12')  # Purine ring
    pyrimidine_pattern = Chem.MolFromSmarts('c1cncnc1')    # Pyrimidine ring

    # Check for nucleobase
    has_purine = mol.HasSubstructMatch(purine_pattern)
    has_pyrimidine = mol.HasSubstructMatch(pyrimidine_pattern)
    if not (has_purine or has_pyrimidine):
        return False, "No nucleobase found"

    # Define SMARTS pattern for pentose sugar (ribose or deoxyribose)
    sugar_pattern = Chem.MolFromSmarts('[C@H]1([O])[C@@H]([O])[C@H]([O])[C@@H](CO)O1')  # Furanose ring
    has_sugar = mol.HasSubstructMatch(sugar_pattern)
    if not has_sugar:
        return False, "No pentose sugar found"

    # Define SMARTS pattern for nucleoside (nucleobase-sugar linkage)
    nucleoside_pattern = Chem.MolFromSmarts('[nH]1[c,nH][c,n][c,n][c,n][c,n]1[C@H]2O[C@@H]([O])[C@H]([O])[C@@H](CO)O2')
    has_nucleoside = mol.HasSubstructMatch(nucleoside_pattern)
    if not has_nucleoside:
        return False, "No nucleoside linkage found"

    # Define SMARTS pattern for phosphate group linked to 3' or 5' hydroxyl group
    phosphate_pattern = Chem.MolFromSmarts('O[P](=O)(O)[O]C[C@H]1O[C@@H]([O])[C@H]([O])[C@@H](O1)CO')  # Phosphate ester
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)
    if not has_phosphate:
        # Check for phosphate at 3' position
        phosphate_3prime_pattern = Chem.MolFromSmarts('O[P](=O)(O)[O][C@H]1O[C@H](CO)[C@@H]([O])[C@H]([O])C1')
        has_phosphate_3prime = mol.HasSubstructMatch(phosphate_3prime_pattern)
        if not has_phosphate_3prime:
            return False, "No phosphate group attached to the sugar's 3' or 5' hydroxyl group"

    return True, "Molecule is a nucleotide with nucleobase, sugar, and phosphate group"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33504',
                              'name': 'nucleotide',
                              'definition': 'A nucleoside phosphate resulting from the condensation of the 3\' or 5\' hydroxy group of a nucleoside with phosphoric acid.',
                              'parents': ['CHEBI:26161']},
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
    'stdout': None}