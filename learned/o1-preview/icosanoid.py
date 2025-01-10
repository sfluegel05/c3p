"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: CHEBI:25099 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    An icosanoid is any signaling molecule arising from oxidation of the three C20 essential fatty acids
    (EFAs): icosapentaenoic acid (EPA), arachidonic acid (AA), and dihomo-gamma-linolenic acid (DGLA).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Since icosanoids are derived from C20 fatty acids, but may have modifications,
    # accept molecules with carbon counts between 18 and 40
    if c_count < 18 or c_count > 40:
        return False, f"Molecule has {c_count} carbon atoms; expected approximately 20 for an icosanoid"

    # Count the number of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Molecule has insufficient oxygen atoms; expected oxidation products"

    # Check for presence of common functional groups
    functional_groups = []

    # Hydroxyl group (OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if hydroxyl_matches:
        functional_groups.append("hydroxyl")

    # Ketone group (>C=O)
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[#6]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    if ketone_matches:
        functional_groups.append("ketone")

    # Carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OX1H1]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if carboxylic_acid_matches:
        functional_groups.append("carboxylic acid")

    # Ester group (-C(=O)O-C)
    ester_pattern = Chem.MolFromSmarts("C(=O)O[CX4]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if ester_matches:
        functional_groups.append("ester")

    # Amide group (-C(=O)N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if amide_matches:
        functional_groups.append("amide")

    # Sulfide group (-S-)
    sulfide_pattern = Chem.MolFromSmarts("S")
    sulfide_matches = mol.GetSubstructMatches(sulfide_pattern)
    if sulfide_matches:
        functional_groups.append("sulfide")

    # Epoxide group (three-membered ring with an oxygen)
    epoxide_pattern = Chem.MolFromSmarts("[C]1[O][C]1")
    epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)
    if epoxide_matches:
        functional_groups.append("epoxide")

    # Peroxide group (-O-O-)
    peroxide_pattern = Chem.MolFromSmarts("[OX2][OX2]")
    peroxide_matches = mol.GetSubstructMatches(peroxide_pattern)
    if peroxide_matches:
        functional_groups.append("peroxide")

    if not functional_groups:
        return False, "Molecule lacks characteristic functional groups of icosanoids"

    # Check for long hydrocarbon chain (saturated or unsaturated)
    # Look for chain of at least 15 carbons
    chain_pattern = Chem.MolFromSmarts("[C][C][C][C][C][C][C][C][C][C][C][C][C][C][C]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Molecule does not have a long hydrocarbon chain"

    # Check for multiple unsaturations (double bonds)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) < 2:
        return False, "Molecule has insufficient double bonds; expected multiple unsaturations"

    return True, f"Molecule has characteristics of an icosanoid (carbon count: {c_count}, functional groups: {', '.join(functional_groups)})"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:25099',
        'name': 'icosanoid',
        'definition': 'Any member of the group of signalling molecules arising from oxidation of the three C20 essential fatty acids (EFAs) icosapentaenoic acid (EPA), arachidonic acid (AA) and dihomo-gamma-linolenic acid (DGLA).',
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
    'attempt': 2,
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