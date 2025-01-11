"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
"""
Classifies: CHEBI:83920 non-proteinogenic amino acid
"""
from rdkit import Chem

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    A non-proteinogenic amino acid is any amino acid that is not naturally encoded in the genetic code of any organism.
    This function will check if the molecule is an amino acid and is not one of the 20 standard amino acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a non-proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for alpha amino acid (N-C-C(=O)-O)
    amino_acid_pattern = Chem.MolFromSmarts('[NX3;H2,H1][CX4][CX3](=O)[O;H1,H0-]')
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Molecule is not an alpha amino acid"

    # Remove stereochemistry from input molecule
    mol_no_stereo = Chem.Mol(mol)
    Chem.RemoveStereochemistry(mol_no_stereo)
    input_can_smiles = Chem.MolToSmiles(mol_no_stereo, isomericSmiles=False)

    # Set of canonical SMILES of standard amino acids without stereochemistry
    standard_amino_acids_smiles = [
        'NCC(=O)O',  # Glycine
        'NCC(C)C(=O)O',  # Alanine
        'NCC(CC)C(=O)O',  # Valine
        'NCC(CC(C)C)C(=O)O',  # Leucine
        'NCC(C(C)C)C(=O)O',  # Isoleucine
        'NCC(CO)C(=O)O',  # Serine
        'NCC(CS)C(=O)O',  # Cysteine
        'NCC(CC(=O)O)C(=O)O',  # Aspartic acid
        'NCC(CC(=O)N)C(=O)O',  # Asparagine
        'NCC(CCC(=O)O)C(=O)O',  # Glutamic acid
        'NCC(CCC(=O)N)C(=O)O',  # Glutamine
        'NCC(CCCNC(N)=N)C(=O)O',  # Arginine
        'NCC(CCCCN)C(=O)O',  # Lysine
        'NCC(CC1=CC=CC=C1)C(=O)O',  # Phenylalanine
        'N1CCCN1C(=O)O',  # Proline
        'NCC(C(O)C)C(=O)O',  # Threonine
        'NCC(CCSC)C(=O)O',  # Methionine
        'NCC(CC1=CNC=N1)C(=O)O',  # Histidine
        'NCC(CC1=CC=CC=C1O)C(=O)O',  # Tyrosine
        'NCC(CC1=CNC2=CC=CC=C12)C(=O)O',  # Tryptophan
    ]

    # Generate canonical SMILES without stereochemistry for standard amino acids
    standard_amino_acids_canonical_smiles = set()
    for std_smiles in standard_amino_acids_smiles:
        std_mol = Chem.MolFromSmiles(std_smiles)
        if std_mol is None:
            continue
        Chem.RemoveStereochemistry(std_mol)
        std_can_smiles = Chem.MolToSmiles(std_mol, isomericSmiles=False)
        standard_amino_acids_canonical_smiles.add(std_can_smiles)

    # Check if input molecule matches any standard amino acid
    if input_can_smiles in standard_amino_acids_canonical_smiles:
        return False, "Molecule is a standard amino acid"

    return True, "Molecule is a non-proteinogenic amino acid"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:83920',
        'name': 'non-proteinogenic amino acid',
        'definition': 'Any amino-acid that is not naturally encoded in the genetic code of any organism.',
        'parents': ['CHEBI:33675']  # parent class is amino acid
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
    'precision': 0.974,
    'recall': 0.867,
    'f1': 0.917,
    'accuracy': 0.9998
}