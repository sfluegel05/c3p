"""
Classifies: CHEBI:64985 bioconjugate
"""
"""
Classifies: CHEBI:XXXXX bioconjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate is defined as a molecular entity consisting of at least 2 biological molecules covalently linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns for biological molecules
    amino_acid_pattern = Chem.MolFromSmarts("[NH2,NH][CX4H][CX3](=O)[OH]")  # Amino acid pattern
    peptide_pattern = Chem.MolFromSmarts("[NH][CX4H][CX3](=O)[NH]")  # Peptide bond pattern
    nucleotide_pattern = Chem.MolFromSmarts("[NX3][CX4][CX4][OX2][CX3](=O)")  # Nucleotide pattern
    carbohydrate_pattern = Chem.MolFromSmarts("[CX4H][OX2H][CX4H][OX2H][CX4H][OX2H]")  # Carbohydrate pattern
    lipid_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]")  # Lipid pattern (long carbon chain)

    # Count matches for each biological molecule pattern
    amino_acid_matches = len(mol.GetSubstructMatches(amino_acid_pattern))
    peptide_matches = len(mol.GetSubstructMatches(peptide_pattern))
    nucleotide_matches = len(mol.GetSubstructMatches(nucleotide_pattern))
    carbohydrate_matches = len(mol.GetSubstructMatches(carbohydrate_pattern))
    lipid_matches = len(mol.GetSubstructMatches(lipid_pattern))

    # Check if at least two biological molecules are present
    biological_molecules_count = 0
    if amino_acid_matches > 0:
        biological_molecules_count += 1
    if peptide_matches > 0:
        biological_molecules_count += 1
    if nucleotide_matches > 0:
        biological_molecules_count += 1
    if carbohydrate_matches > 0:
        biological_molecules_count += 1
    if lipid_matches > 0:
        biological_molecules_count += 1

    if biological_molecules_count < 2:
        return False, f"Only {biological_molecules_count} biological molecule(s) found, need at least 2"

    # Check for covalent linkage (e.g., ester, amide, thioether bonds)
    covalent_linkage_pattern = Chem.MolFromSmarts("[CX4][O,N,S][CX4]")  # Covalent linkage pattern
    covalent_linkage_matches = len(mol.GetSubstructMatches(covalent_linkage_pattern))
    if covalent_linkage_matches == 0:
        return False, "No covalent linkage found between biological molecules"

    return True, "Contains at least 2 biological molecules covalently linked together"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:XXXXX',
                          'name': 'bioconjugate',
                          'definition': 'A molecular entity consisting of at least 2 biological molecules covalently linked together.',
                          'parents': ['CHEBI:XXXXX', 'CHEBI:XXXXX']},
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