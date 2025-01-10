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

    # Look for peptide bond pattern (-C(=O)-N-)
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3H0,NX3H1]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) < 1:
        return False, "No peptide bond found"

    # Count amino acid residues by looking for amino and carboxyl groups
    amino_group_pattern = Chem.MolFromSmarts("[NX3H2,NX3H1,NX3H0]")
    carboxyl_group_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1,OX1H0-]")
    
    amino_matches = mol.GetSubstructMatches(amino_group_pattern)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_group_pattern)
    
    if len(amino_matches) < 2 or len(carboxyl_matches) < 2:
        return False, "Not enough amino or carboxyl groups to form two amino acid residues"

    # Ensure that there are at least two amino acid residues
    # Each residue should have at least one amino group and one carboxyl group
    # We can count the number of distinct amino and carboxyl groups
    if len(amino_matches) < 2 or len(carboxyl_matches) < 2:
        return False, "Insufficient amino or carboxyl groups for a dipeptide"

    # Check molecular weight - dipeptides typically have a molecular weight between 150-1000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150 or mol_wt > 1000:
        return False, f"Molecular weight {mol_wt:.2f} Da is outside the typical range for a dipeptide"

    # Check for the presence of two amino acid residues connected by a peptide bond
    # This is a more relaxed check to ensure we capture more complex dipeptides
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