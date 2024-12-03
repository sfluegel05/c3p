"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide (contains four amino-acid residues connected by peptide linkages).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define peptide bond pattern (C(=O)N)
    peptide_bond_pattern = Chem.MolFromSmarts('C(=O)N')

    # Find all peptide bonds
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)

    # Check if there are exactly 3 peptide bonds (4 amino acids connected by 3 peptide bonds)
    if len(peptide_bonds) != 3:
        return False, f"Expected 3 peptide bonds, found {len(peptide_bonds)}"

    # Define amino acid residue pattern (N-C(alpha)-C(=O) pattern)
    amino_acid_pattern = Chem.MolFromSmarts('[NX3][CX4][CX3](=O)')
    amino_acid_residues = mol.GetSubstructMatches(amino_acid_pattern)

    # Check if there are exactly 4 amino acid residues
    if len(amino_acid_residues) != 4:
        return False, f"Expected 4 amino acid residues, found {len(amino_acid_residues)}"

    return True, "Valid tetrapeptide"

# Example usage
smiles = "CC(C)C[C@H](NC(=O)[C@H](C)N)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CO)C(O)=O"
print(is_tetrapeptide(smiles))  # Should return (True, "Valid tetrapeptide")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:48030',
                          'name': 'tetrapeptide',
                          'definition': 'Any molecule that contains four '
                                        'amino-acid residues connected by '
                                        'peptide linkages.',
                          'parents': ['CHEBI:25676']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Valid tetrapeptide')\n",
    'num_true_positives': 7,
    'num_false_positives': 0,
    'num_true_negatives': 13,
    'num_false_negatives': 6,
    'precision': 1.0,
    'recall': 0.5384615384615384,
    'f1': 0.7000000000000001,
    'accuracy': None}