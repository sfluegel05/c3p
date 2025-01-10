"""
Classifies: CHEBI:33838 nucleoside
"""
from rdkit import Chem

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside consists of a nitrogenous base linked to a sugar moiety via a β-glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for identifying purine and pyrimidine bases
    purine_pattern = Chem.MolFromSmarts('c1ncnc2[nH]cnc12')  # Adjusted to catch common variations
    pyrimidine_pattern = Chem.MolFromSmarts('c1ccncn1')      # Generalized pyrimidine structure
    if not (mol.HasSubstructMatch(purine_pattern) or mol.HasSubstructMatch(pyrimidine_pattern)):
        return False, "No recognizable nucleobase pattern found"
    
    # Generalized SMARTS pattern for ribose or deoxyribose units
    sugar_pattern = Chem.MolFromSmarts('O[C@@H]1[C@H](O)[C@@H](O)[C@H](CO)O1')  # Capture flexibility in sugar structure
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No recognizable ribose or deoxyribose pattern found"
    
    # Validate glycosidic bond pattern
    glycosidic_bond_pattern = Chem.MolFromSmarts('CO(c1[nH]c(cn1))')  # Match `base-O-sugar` linkage
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No proper glycosidic bond between sugar and base"

    return True, "Contains nucleobase and sugar moiety with appropriate glycosidic bond"


# Metadata and configuration
__metadata__ = {
    'chemical_class': {
        'id': 'None',
        'name': 'nucleoside',
        'definition': 'Nucleosides are compounds formed by linking a nucleobase to a sugar moiety through a glycosidic bond.',
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5
    }
}