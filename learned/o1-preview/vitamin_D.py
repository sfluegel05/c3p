"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: vitamin D
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D based on its SMILES string.
    Vitamin D refers to a group of fat-soluble hydroxy seco-steroids that exhibit biological activity against vitamin D deficiency.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is vitamin D, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define vitamin D core skeleton SMARTS pattern
    # This pattern represents the secosteroid skeleton common to vitamin D molecules
    vitamin_D_skeleton_smarts = 'C1=CC[C@H]2C1CCC3C2=C\C=C\C4=CCCCC34'
    skeleton_pattern = Chem.MolFromSmarts(vitamin_D_skeleton_smarts)
    if not mol.HasSubstructMatch(skeleton_pattern):
        return False, "No vitamin D secosteroid skeleton found"

    # Check for at least one hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts('[OX2H]')
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) < 1:
        return False, "No hydroxy groups found"

    # Check that molecule is lipophilic (logP value)
    logP = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    if logP < 3:
        return False, "LogP value too low, molecule may not be fat-soluble"

    # Check molecular weight - vitamin D molecules typically ~380-420 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 450:
        return False, "Molecular weight not in typical range for vitamin D"

    return True, "Molecule matches vitamin D core skeleton with hydroxy groups"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'vitamin D',
        'definition': 'Any member of a group of fat-soluble hydroxy seco-steroids that exhibit biological activity against vitamin D deficiency. Vitamin D can be obtained from sun exposure, food and supplements and is biologically inactive and converted into the biologically active calcitriol via double hydroxylation in the body.',
        'parents': []
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
    'stdout': None
}