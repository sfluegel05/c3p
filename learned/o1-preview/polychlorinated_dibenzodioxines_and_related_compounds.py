"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
To determine if a molecule is a polychlorinated dibenzodioxin, dibenzofuran, biphenyl, or a related compound, I will:

1. Parse the SMILES string into an RDKit molecule.
2. Check for the presence of the dibenzodioxin, dibenzofuran, or biphenyl core structures using substructure matching.
3. Count the number of chlorine (Cl) and bromine (Br) atoms in the molecule.
4. Verify that the molecule is polychlorinated or polybrominated by checking that there are at least 2 halogen atoms.
5. Return `True` if all conditions are met, along with a reason for the classification.

Here is the code:

"""
Classifies: polychlorinated dibenzodioxines and related compounds

"""
from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule is a polychlorinated dibenzodioxin, dibenzofuran, biphenyl, or a related compound based on its SMILES string.

    The compound must contain one of the core structures (dibenzodioxin, dibenzofuran, or biphenyl) and must be polychlorinated or polybrominated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polychlorinated dibenzodioxin or related compound, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define core structures
    biphenyl_pattern = Chem.MolFromSmarts('c1ccccc1-c2ccccc2')  # Biphenyl core
    dibenzodioxin_pattern = Chem.MolFromSmiles('O1c2ccccc2Oc3ccccc13')  # Dibenzodioxin core
    dibenzofuran_pattern = Chem.MolFromSmiles('c1cc2oc3ccccc3cc2cc1')  # Dibenzofuran core

    # Check for core structures
    has_biphenyl = mol.HasSubstructMatch(biphenyl_pattern)
    has_dibenzodioxin = mol.HasSubstructMatch(dibenzodioxin_pattern)
    has_dibenzofuran = mol.HasSubstructMatch(dibenzofuran_pattern)

    if not (has_biphenyl or has_dibenzodioxin or has_dibenzofuran):
        return False, "Does not contain dibenzodioxin, dibenzofuran, or biphenyl core structure"

    # Count Cl and Br atoms
    num_Cl = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17)
    num_Br = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 35)

    num_halogen = num_Cl + num_Br

    if num_halogen < 2:
        return False, f"Contains {num_halogen} halogen atoms, needs at least 2 to be polychlorinated or polybrominated"

    return True, "Contains polychlorinated/brominated dibenzodioxin, dibenzofuran, or biphenyl core"


__metadata__ = {'chemical_class': {
    'id': None,
    'name': 'polychlorinated dibenzodioxines and related compounds',
    'definition': 'Organochlorine compounds that are polychlorinated dibenzodioxines and structurally related entities that are persistent organic pollutants. These include polychlorinated dibenzofurans as well as polychlorinated and polybrominated biphenyls  They vary widely in their toxicity, but their toxic mode of action is through the aryl hydrocarbon receptor.',
    'parents': []
    },
    'config': {   'llm_model_name': 'YourModelNameHere',
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

"""