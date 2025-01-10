"""
Classifies: CHEBI:33856 aromatic amino acid
"""
"""
Classifies: CHEBI:33857 aromatic amino acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid is defined as an amino acid whose structure includes an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for amino acid pattern: at least one amino group (NH2 or NH) and one carboxyl group (COOH)
    amino_pattern = Chem.MolFromSmarts("[NH2,NH]")
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # Check for aromatic rings
    aromatic_rings = mol.GetSubstructMatches(Chem.MolFromSmarts("[a]"))
    if not aromatic_rings:
        return False, "No aromatic ring found"

    # Ensure the aromatic ring is part of the molecule (not just a fragment)
    # and is connected to the amino acid backbone
    for ring in aromatic_rings:
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == "C" and neighbor.GetDegree() > 1:
                    return True, "Contains amino acid functional groups and an aromatic ring"

    return False, "Aromatic ring not connected to amino acid backbone"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:33857',
        'name': 'aromatic amino acid',
        'definition': 'An amino acid whose structure includes an aromatic ring.',
        'parents': ['CHEBI:33704', 'CHEBI:33857']
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