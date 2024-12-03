"""
Classifies: CHEBI:83822 non-proteinogenic L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

proteinogenic_amino_acids = {
    "Glycine": "C(C(=O)O)N",
    "Alanine": "CC(C(=O)O)N",
    "Valine": "CC(C)C(C(=O)O)N",
    "Leucine": "CC(C)CC(C(=O)O)N",
    "Isoleucine": "CC(C)C(C(=O)O)N",
    "Proline": "C1CC(C1)C(=O)O",
    "Phenylalanine": "C1=CC=C(C=C1)CC(C(=O)O)N",
    "Tyrosine": "C1=CC(=CC=C1CC(C(=O)O)N)O",
    "Tryptophan": "C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N",
    "Serine": "C(C(C(=O)O)N)O",
    "Threonine": "C(C(C(C(=O)O)N)O)C",
    "Cysteine": "C(C(C(=O)O)N)S",
    "Methionine": "CC(C)CC(C(=O)O)N",
    "Asparagine": "C(C(C(=O)O)N)C(=O)N",
    "Glutamine": "C(CC(C(=O)O)N)C(=O)N",
    "Aspartic acid": "C(C(C(=O)O)N)C(=O)O",
    "Glutamic acid": "C(CC(C(=O)O)N)C(=O)O",
    "Lysine": "C(CCN)CC(C(=O)O)N",
    "Arginine": "C(C(C(=O)O)N)NC(=N)N",
    "Histidine": "C1=CC(=CN1)CC(C(=O)O)N",
    "Selenocysteine": "C(C(C(=O)O)N)Se",
    "Pyrrolysine": "C1CC(C1)C(C(=O)O)N",
}

def is_non_proteinogenic_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic L-alpha-amino acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a non-proteinogenic L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if it has the general structure of an alpha-amino acid
    alpha_amino_acid = Chem.MolFromSmarts("N[C@@H](C(=O)O)C")
    if not mol.HasSubstructMatch(alpha_amino_acid):
        return False, "Does not match the general structure of an L-alpha-amino acid"

    # Check if it is one of the 23 proteinogenic amino acids
    for name, prot_smiles in proteinogenic_amino_acids.items():
        prot_mol = Chem.MolFromSmiles(prot_smiles)
        if mol.HasSubstructMatch(prot_mol):
            return False, f"Matches the structure of the proteinogenic amino acid: {name}"

    return True, "Non-proteinogenic L-alpha-amino acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:83822',
                          'name': 'non-proteinogenic L-alpha-amino acid',
                          'definition': 'Any L-alpha-amino acid which is not a '
                                        'member of the group of 23 '
                                        'proteinogenic amino acids.',
                          'parents': ['CHEBI:15705', 'CHEBI:83925']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 2-3: malformed \\N character escape (<string>, line 1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}