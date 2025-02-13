"""
Classifies: CHEBI:48030 tetrapeptide
"""
"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide is a molecule that contains four amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the peptide bond pattern (-C(=O)-N-)
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3H0,H1]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # A tetrapeptide should have at least 3 peptide bonds (connecting 4 amino acids)
    if len(peptide_bond_matches) < 3:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, need at least 3"

    # Count the number of amino acid residues using a more flexible pattern
    # Each amino acid residue should have at least a carboxyl group (-C(=O)-) and an amino group (-NH-)
    amino_acid_pattern = Chem.MolFromSmarts("[CX4][CX3](=[OX1])[OX2H0,H1].[NX3H0,H1,H2]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    
    # A tetrapeptide should have at least 4 amino acid residues
    if len(amino_acid_matches) < 4:
        return False, f"Found {len(amino_acid_matches)} amino acid-like groups, need at least 4"

    # Check for the presence of at least one terminal carboxyl group and one terminal amino group
    terminal_carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    terminal_amino_pattern = Chem.MolFromSmarts("[NX3H2]")
    
    terminal_carboxyl_matches = mol.GetSubstructMatches(terminal_carboxyl_pattern)
    terminal_amino_matches = mol.GetSubstructMatches(terminal_amino_pattern)
    
    if len(terminal_carboxyl_matches) < 1 and len(terminal_amino_matches) < 1:
        return False, "Missing both terminal carboxyl and amino groups"

    # Check molecular weight - tetrapeptides typically have a molecular weight between 300-1000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1000:
        return False, f"Molecular weight {mol_wt:.2f} Da is outside the typical range for tetrapeptides"

    # Check for linearity of the peptide chain
    # Count the number of rotatable bonds - should be at least 3 for a tetrapeptide
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Too few rotatable bonds for a tetrapeptide"

    return True, "Contains four amino acid residues connected by at least three peptide bonds"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:48030',
                          'name': 'tetrapeptide',
                          'definition': 'Any molecule that contains four amino-acid residues connected by peptide linkages.',
                          'parents': ['CHEBI:16670', 'CHEBI:26561']},
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