"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: CHEBI:25676 oligopeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is a peptide containing a relatively small number of amino acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligopeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for peptide bonds (amide bonds, -C(=O)-N-)
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) == 0:
        return False, "No peptide bonds found"

    # Count the number of peptide bonds to estimate the number of amino acids
    # Each peptide bond connects two amino acids, so the number of amino acids is roughly (number of peptide bonds + 1)
    num_amino_acids = len(peptide_bond_matches) + 1

    # Oligopeptides typically have 2 to 20 amino acids
    if num_amino_acids < 2:
        return False, f"Too few amino acids: {num_amino_acids}"
    if num_amino_acids > 20:
        return False, f"Too many amino acids: {num_amino_acids}"

    # Check for typical amino acid side chains (e.g., -CH3, -CH2-, -COOH, -NH2, etc.)
    side_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    side_chain_matches = mol.GetSubstructMatches(side_chain_pattern)
    if len(side_chain_matches) < num_amino_acids:
        return False, f"Missing side chains for {num_amino_acids} amino acids"

    # Check for the presence of both amino and carboxyl groups (typical of amino acids)
    amino_group_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")
    carboxyl_group_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    amino_matches = mol.GetSubstructMatches(amino_group_pattern)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_group_pattern)
    if len(amino_matches) == 0 or len(carboxyl_matches) == 0:
        return False, "Missing amino or carboxyl groups"

    return True, f"Contains {num_amino_acids} amino acids linked by peptide bonds"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:25676',
        'name': 'oligopeptide',
        'definition': 'A peptide containing a relatively small number of amino acids.',
        'parents': ['CHEBI:16670', 'CHEBI:26561']
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
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199
}