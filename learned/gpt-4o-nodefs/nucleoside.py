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
    purine_pattern = Chem.MolFromSmarts('n1cnc2c1ncnc2')  # Purine base pattern
    pyrimidine_pattern = Chem.MolFromSmarts('n1c([nH])cnc1=O')  # Simplified pyrimidine pattern
    if not (mol.HasSubstructMatch(purine_pattern) or mol.HasSubstructMatch(pyrimidine_pattern)):
        return False, "No recognizable nucleobase pattern found"
    
    # Generalized SMARTS pattern for ribose or deoxyribose units
    sugar_pattern = Chem.MolFromSmarts('O[C@@H]1[C@H](O)[C@@H](O)[C@H](CO)O1')  # Capture flexibility in sugar structure
    alternate_sugar_pattern = Chem.MolFromSmarts('O[C@H]1[C@@H](O)[C@@H](O)[C@H](CO)O1')  # Another orientation of ribose
    if not (mol.HasSubstructMatch(sugar_pattern) or mol.HasSubstructMatch(alternate_sugar_pattern)):
        return False, "No recognizable ribose or deoxyribose pattern found"
    
    # Validate presence of a glycosidic bond between a base and a sugar
    # We assume the link involves a nitrogen in the base connecting to the anomeric carbon of the sugar
    glycosidic_bond_found = False
    for bond in mol.GetBonds():
        if (bond.GetBeginAtom().GetAtomicNum() == 8 and bond.GetEndAtom().GetAtomicNum() == 6) or \
           (bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 8):
            if (mol.HasSubstructMatch(purine_pattern) or mol.HasSubstructMatch(pyrimidine_pattern)) and \
               (mol.HasSubstructMatch(sugar_pattern) or mol.HasSubstructMatch(alternate_sugar_pattern)):
                glycosidic_bond_found = True
                break
    
    if not glycosidic_bond_found:
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