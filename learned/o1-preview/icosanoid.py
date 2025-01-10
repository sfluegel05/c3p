"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: CHEBI:25099 icosanoid
"""
from rdkit import Chem

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

    # Since icosanoids are derived from C20 fatty acids, accept molecules with carbon counts between 18 and 25
    if c_count < 18 or c_count > 25:
        return False, f"Molecule has {c_count} carbon atoms; expected approximately 20 for an icosanoid"

    # Count the number of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Molecule has insufficient oxygen atoms; expected oxidation products"

    # Check for presence of common functional groups
    functional_groups = []

    # Hydroxyl group (OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if mol.HasSubstructMatch(hydroxyl_pattern):
        functional_groups.append("hydroxyl")

    # Ketone group (>C=O)
    ketone_pattern = Chem.MolFromSmarts("C(=O)[C]")
    if mol.HasSubstructMatch(ketone_pattern):
        functional_groups.append("ketone")

    # Aldehyde group (-CHO)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    if mol.HasSubstructMatch(aldehyde_pattern):
        functional_groups.append("aldehyde")

    # Carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O[H]")
    if mol.HasSubstructMatch(carboxylic_acid_pattern):
        functional_groups.append("carboxylic acid")

    # Epoxide group (three-membered ring with an oxygen)
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")
    if mol.HasSubstructMatch(epoxide_pattern):
        functional_groups.append("epoxide")

    # Peroxide group (-O-O-)
    peroxide_pattern = Chem.MolFromSmarts("O[O]")
    if mol.HasSubstructMatch(peroxide_pattern):
        functional_groups.append("peroxide")

    # Ether group (-O-)
    ether_pattern = Chem.MolFromSmarts("C-O-C")
    if mol.HasSubstructMatch(ether_pattern):
        functional_groups.append("ether")

    # Lactone group (cyclic ester)
    lactone_pattern = Chem.MolFromSmarts("C1OC(=O)C1")
    if mol.HasSubstructMatch(lactone_pattern):
        functional_groups.append("lactone")

    # If no functional groups found, return False
    if not functional_groups:
        return False, "Molecule lacks characteristic functional groups of icosanoids"

    # Since icosanoids are oxidation products of fatty acids,
    # they typically have at least one oxygen-containing functional group
    # and a hydrophobic hydrocarbon chain.

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
    'attempt': 3,
    'success': True,
    'best': True
}