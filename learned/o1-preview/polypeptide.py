"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: polypeptide
"""
from rdkit import Chem

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

    # Define peptide bond pattern (amide bond between amino acids)
    peptide_bond_smarts = "[NX3][CX3](=O)"
    peptide_bond_pattern = Chem.MolFromSmarts(peptide_bond_smarts)
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bond_matches)
    
    # Estimate number of amino acid residues
    num_residues = num_peptide_bonds + 1  # Number of residues is number of peptide bonds plus one
    
    if num_residues >= 10:
        return True, f"Contains {num_residues} amino acid residues"
    else:
        return False, f"Contains only {num_residues} amino acid residues, which is less than 10"

__metadata__ = {  
    'chemical_class': {  
        'id': None,
        'name': 'polypeptide',
        'definition': 'A peptide containing ten or more amino acid residues.',
        'parents': None
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
    # Additional metadata fields can be added here
}