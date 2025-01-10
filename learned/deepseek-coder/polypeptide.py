"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: CHEBI:16670 polypeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide based on its SMILES string.
    A polypeptide is a peptide containing ten or more amino acid residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for peptide bonds (amide bonds) using SMARTS pattern
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3H1]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # Each peptide bond corresponds to one amino acid residue
    num_residues = len(peptide_bond_matches) + 1  # N-1 bonds for N residues
    
    if num_residues < 10:
        return False, f"Found {num_residues} amino acid residues, need at least 10"
    
    # Check for typical polypeptide features (e.g., multiple amide bonds, large molecular weight)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:  # Polypeptides are typically large molecules
        return False, "Molecular weight too low for polypeptide"
    
    # Count nitrogen and oxygen atoms (polypeptides have many due to amide bonds)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_count < 10 or o_count < 10:
        return False, "Too few nitrogen or oxygen atoms for polypeptide"
    
    return True, f"Contains {num_residues} amino acid residues, qualifies as a polypeptide"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:16670',
        'name': 'polypeptide',
        'definition': 'A peptide containing ten or more amino acid residues.',
        'parents': ['CHEBI:16670', 'CHEBI:16670']
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