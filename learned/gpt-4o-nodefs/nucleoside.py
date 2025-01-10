"""
Classifies: CHEBI:33838 nucleoside
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside typically consists of a nitrogenous base linked to a sugar moiety via a β-glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for identifying nucleobases
    purine_pattern = Chem.MolFromSmarts('n1(cnc2c1ncnc2)n')
    pyrimidine_pattern = Chem.MolFromSmarts('c1c[nH]c(=O)[nH]c1=O')

    # Check for presence of nucleobase
    has_purine = mol.HasSubstructMatch(purine_pattern)
    has_pyrimidine = mol.HasSubstructMatch(pyrimidine_pattern)
    if not (has_purine or has_pyrimidine):
        return False, "No nucleobase found"

    # SMARTS pattern for ribose or deoxyribose
    ribose_pattern = Chem.MolFromSmarts('O[C@@H]1[C@H](O)[C@@H](O)[C@H](CO)O1')
    deoxyribose_pattern = Chem.MolFromSmarts('O[C@@H]1[C@H](O)[C@H](O)[C@H](CO)O1')

    # Check for sugar moiety
    has_ribose = mol.HasSubstructMatch(ribose_pattern)
    has_deoxyribose = mol.HasSubstructMatch(deoxyribose_pattern)
    if not (has_ribose or has_deoxyribose):
        return False, "No ribose or deoxyribose found"
    
    # Check for glycosidic bond connectivity
    # Example glycosidic bond pattern check (simplified for illustration)
    glycosidic_bond_pattern = Chem.MolFromSmarts('n-C-O')
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond between sugar and base"

    return True, "Contains nucleobase and sugar moiety with appropriate glycosidic bond"


# Metadata and configuration
__metadata__ = {
    'chemical_class': {
        'id': 'None',  # Replace with actual ID if available
        'name': 'nucleoside',
        'definition': 'Nucleosides are compounds formed by linking a nucleobase to a sugar moiety through a glycosidic bond.',
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5
    }
}