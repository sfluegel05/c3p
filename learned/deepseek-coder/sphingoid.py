"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: CHEBI:37684 sphingoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    A sphingoid is characterized by a long hydrocarbon chain with a hydroxyl group,
    an amino group, and often a double bond or additional hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a long hydrocarbon chain (at least 10 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Molecule has too few carbons to be a sphingoid"

    # Check for the presence of a hydroxyl group (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Check for the presence of an amino group (-NH2, -NH3+, or -NH-)
    amino_pattern = Chem.MolFromSmarts("[NX3H2,NX4H3+,NX3H]")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"

    # Check for the presence of a 2-amino-1,3-diol or 2-amino-1,3,4-triol backbone
    backbone_pattern = Chem.MolFromSmarts("[CX4][CX4]([OH])[CX4]([NH2,NH3+,NH])[CX4]([OH])")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No 2-amino-1,3-diol or 2-amino-1,3,4-triol backbone found"

    # Check for the presence of a double bond (optional, as some sphingoids are saturated)
    double_bond_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    has_double_bond = mol.HasSubstructMatch(double_bond_pattern)

    # Check for additional hydroxyl groups (optional, as some sphingoids have more than one)
    additional_hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(additional_hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "No additional hydroxyl groups found"

    return True, "Contains a long hydrocarbon chain with a hydroxyl group, an amino group, and a 2-amino-1,3-diol or 2-amino-1,3,4-triol backbone"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37684',
                          'name': 'sphingoid',
                          'definition': 'Sphinganine, its homologs and stereoisomers, and the hydroxy and unsaturated derivatives of these compounds.',
                          'parents': ['CHEBI:37684']},
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
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}