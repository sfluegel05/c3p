"""
Classifies: CHEBI:33916 aldopentose
"""
from rdkit import Chem

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.

    An aldopentose is defined as a pentose (5-carbon sugar) with a potential aldehyde
    group at one end when in linear form. Often exists in cyclic form (furanose or 
    pyranose) within the context of hemiacetal structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aldopentose, False otherwise.
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if there are exactly 5 carbon atoms (pentose sugar)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 5:
        return False, f"Expected 5 carbon atoms, found {c_count}"

    # Check for either linear aldehyde group or cyclic forms typical for aldoses
    # Applying strict pattern for 5-member furanose or 6-member pyranose
    furanose_pattern = Chem.MolFromSmarts("[O][C@@H]1O[C@H]2O[C@H]([C@H]([C@H]1[O])O)[C@H]2")
    pyranose_pattern = Chem.MolFromSmarts("[O][C@H]1O[C@H](O)[C@H](O)[C@H]1")
    
    if mol.HasSubstructMatch(furanose_pattern) or mol.HasSubstructMatch(pyranose_pattern):
        return True, "Contains cyclic furanose or pyranose structure typical for aldopentoses"

    # Check for the presence of an aldehyde group (R-CHO) in linear forms
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Contains a 5-carbon backbone with an aldehyde group at one end"

    return False, "Does not match aldopentose structure requirements"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:26977',
        'name': 'aldopentose',
        'definition': 'A pentose with a (potential) aldehyde group at one end.',
        'parents': ['CHEBI:4705']
    },
    'config': {
        'llm_model_name': 'your-model',
        'f1_threshold': 0.8,
        'max_attempts': 5
    },
    'message': None,
    'attempt': 0,
    'success': True
}