"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
"""
Classifies: nitrohydrocarbon
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon in which one or more hydrogens has been replaced by nitro groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check that molecule contains only C, H, N, O atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [6, 1, 7, 8]:
            return False, f"Contains atom other than C, H, N, O: {atom.GetSymbol()}"

    # Identify nitro groups attached to carbon atoms
    nitro_group = Chem.MolFromSmarts('[N+](=O)[O-]')
    carbon_nitro = Chem.MolFromSmarts('[C][N+](=O)[O-]')
    nitro_matches = mol.GetSubstructMatches(carbon_nitro)
    if not nitro_matches:
        return False, "No nitro groups attached to carbon found"

    # Check for functional groups other than nitro groups
    # For simplicity, we can search for common functional groups
    non_nitro_functional_groups = [
        Chem.MolFromSmarts(smarts) for smarts in [
            '[CX3]=[OX1]',       # Carbonyl group
            '[CX3]-[OX2H]',      # Alcohol
            '[CX4][Cl,Br,I,F]',  # Halogens
            '[C]=[C]',           # Alkenes
            '[C]#[C]',           # Alkynes
            '[CX3]-[OX2]-[CX3]', # Ethers
            '[CX3](=O)[OX2H1]',  # Carboxylic acid
            '[CX3](=O)[NX3H2]',  # Amide
            '[NX3H2]',           # Primary amine
            '[NX3H][CX3]',       # Secondary amine
            '[NX3]([CX3])[CX3]', # Tertiary amine
        ]
    ]
    for fg in non_nitro_functional_groups:
        if mol.HasSubstructMatch(fg):
            return False, "Contains other functional groups besides nitro"

    return True, "Molecule is a hydrocarbon with nitro groups attached to carbon atoms"

__metadata__ = {   'chemical_class': {   'id': 'CUSTOM',
                          'name': 'nitrohydrocarbon',
                          'definition': 'A C-nitro compound that is a hydrocarbon in which one or more of the hydrogens has been replaced by nitro groups.',
                          'parents': []},
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
    'accuracy': None}