"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide is a molecule that contains two amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dipeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns for peptide bond and amino acid residues
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3H0,NX3H1]")
    amino_acid_pattern = Chem.MolFromSmarts("[NX3H2,NX3H1,NX3H0][CX4H][CX3](=[OX1])[OX2H1,OX1H0-]")
    
    # Find peptide bonds and amino acid residues
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    
    # Check for at least one peptide bond
    if len(peptide_bond_matches) < 1:
        return False, "No peptide bond found"
    
    # Check for at least two amino acid residues
    if len(amino_acid_matches) < 2:
        return False, "Not enough amino acid residues found"
    
    # Verify that the peptide bond connects two amino acid residues
    # Get the atoms involved in the peptide bond
    peptide_atoms = set()
    for match in peptide_bond_matches:
        peptide_atoms.update(match)
    
    # Check if the peptide bond connects two amino acid residues
    connected_residues = 0
    for match in amino_acid_matches:
        if any(atom in peptide_atoms for atom in match):
            connected_residues += 1
    
    if connected_residues < 2:
        return False, "Peptide bond does not connect two amino acid residues"
    
    # Check molecular weight - dipeptides typically have a molecular weight between 150-1500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150 or mol_wt > 1500:
        return False, f"Molecular weight {mol_wt:.2f} Da is outside the typical range for a dipeptide"
    
    return True, "Contains two amino acid residues connected by a peptide bond"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:46761',
        'name': 'dipeptide',
        'definition': 'Any molecule that contains two amino-acid residues connected by peptide linkages.',
        'parents': ['CHEBI:46761', 'CHEBI:46761']
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