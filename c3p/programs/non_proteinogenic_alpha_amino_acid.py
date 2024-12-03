"""
Classifies: CHEBI:83925 non-proteinogenic alpha-amino acid
"""
from rdkit import Chem

# List of SMILES strings for the 23 proteinogenic amino acids
proteinogenic_amino_acids_smiles = [
    "C(C(=O)O)N",  # Glycine
    "CC(C(=O)O)N",  # Alanine
    "CC(C)C(C(=O)O)N",  # Valine
    "CCC(C)C(C(=O)O)N",  # Leucine
    "CC(C)CC(C(=O)O)N",  # Isoleucine
    "C(C(C(=O)O)N)O",  # Serine
    "CC(C(C(=O)O)N)O",  # Threonine
    "C(C(C(=O)O)N)S",  # Cysteine
    "CC(C(=O)O)N",  # Methionine
    "C1CC[C@H](C1)C(C(=O)O)N",  # Proline
    "C1=CC=C(C=C1)CC(C(=O)O)N",  # Phenylalanine
    "C1=CC(=CC=C1CC(C(=O)O)N)O",  # Tyrosine
    "C1=CC2=C(C=C1)C(=CN2)CC(C(=O)O)N",  # Tryptophan
    "C(C(C(=O)O)N)C(=O)O",  # Aspartic acid
    "CC(C(=O)O)C(=O)O",  # Glutamic acid
    "C(C(C(=O)O)N)C(=O)N",  # Asparagine
    "CCC(C(=O)O)C(=O)N",  # Glutamine
    "C(CC(C(=O)O)N)C(=O)O",  # Lysine
    "C(CC(C(=O)O)N)CC(=O)O",  # Arginine
    "C(C(C(=O)O)N)S",  # Histidine
    "C1=CC=C(C=C1)C(C(=O)O)N",  # Phenylalanine
    "CC(C(=O)O)N",  # Threonine
    "C(C(=O)O)N",  # Glutamine
]

def is_non_proteinogenic_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is a non-proteinogenic alpha-amino acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a non-proteinogenic alpha-amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is an alpha-amino acid
    alpha_amino_acid = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetTotalDegree() == 3:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and neighbor.GetTotalDegree() == 4:
                    for neighbor2 in neighbor.GetNeighbors():
                        if neighbor2.GetSymbol() == 'C' and neighbor2.GetTotalDegree() == 4:
                            alpha_amino_acid = True
                            break
                    if alpha_amino_acid:
                        break
            if alpha_amino_acid:
                break

    if not alpha_amino_acid:
        return False, "Not an alpha-amino acid"

    # Check if the molecule is one of the 23 proteinogenic amino acids
    for proteinogenic_smiles in proteinogenic_amino_acids_smiles:
        proteinogenic_mol = Chem.MolFromSmiles(proteinogenic_smiles)
        if mol.HasSubstructMatch(proteinogenic_mol):
            return False, "Molecule is a proteinogenic amino acid"

    return True, "Molecule is a non-proteinogenic alpha-amino acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:83925',
                          'name': 'non-proteinogenic alpha-amino acid',
                          'definition': 'Any  alpha-amino acid which is not a '
                                        'member of the group of 23 '
                                        'proteinogenic amino acids.',
                          'parents': ['CHEBI:33704', 'CHEBI:83820']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 1-2: malformed \\N character escape (<string>, line 1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}