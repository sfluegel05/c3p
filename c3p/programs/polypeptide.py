"""
Classifies: CHEBI:15841 polypeptide
"""
from rdkit import Chem

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide (a peptide containing ten or more amino acid residues).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define amino acid residue patterns
    amino_acid_patterns = [
        "N[C@@H](C)C(=O)",  # General amino acid backbone
        "N[C@H](C)C(=O)",   # General amino acid backbone (chiral center)
    ]

    # Convert patterns to RDKit molecules
    amino_acid_mols = [Chem.MolFromSmarts(pattern) for pattern in amino_acid_patterns]

    # Find all matches for amino acid residues
    matches = set()
    for aa_mol in amino_acid_mols:
        matches.update(mol.GetSubstructMatches(aa_mol))

    # Count unique amino acid residues
    num_residues = len(matches)

    if num_residues >= 10:
        return True, f"Polypeptide with {num_residues} amino acid residues"
    else:
        return False, f"Only {num_residues} amino acid residues found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:15841',
                          'name': 'polypeptide',
                          'definition': 'A peptide containing ten or more '
                                        'amino acid residues.',
                          'parents': ['CHEBI:16670', 'CHEBI:33839']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 15,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 31,
    'precision': 1.0,
    'recall': 0.32608695652173914,
    'f1': 0.4918032786885246,
    'accuracy': None}