"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
"""
Classifies: CHEBI:58937 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    An N-acyl-L-alpha-amino acid is defined as any L-alpha-amino acid carrying an N-acyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for L-alpha-amino acid backbone: [C@@H](C(=O)O)N
    amino_acid_pattern = Chem.MolFromSmarts("[C@@H](C(=O)O)N")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No L-alpha-amino acid backbone found"

    # Get the nitrogen atom in the amino acid backbone
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    amino_acid_nitrogen = amino_acid_matches[0][2]  # Nitrogen atom index in the amino acid pattern

    # Check for N-acyl substituent: N-C(=O)
    acyl_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No N-acyl substituent found"

    # Get all nitrogen atoms in the molecule
    nitrogen_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]

    # Check if any nitrogen atom is part of an N-acyl group and is connected to the amino acid structure
    for nitrogen in nitrogen_atoms:
        # Check if the nitrogen is part of an N-acyl group
        for neighbor in mol.GetAtomWithIdx(nitrogen).GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE for bond in neighbor.GetBonds()):
                # Check if the nitrogen is part of the amino acid structure (backbone or side chain)
                if nitrogen == amino_acid_nitrogen or mol.GetAtomWithIdx(nitrogen).GetDegree() > 1:
                    return True, "Contains L-alpha-amino acid backbone with N-acyl substituent"

    return False, "N-acyl group is not attached to the amino group of the amino acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:58937',
                          'name': 'N-acyl-L-alpha-amino acid',
                          'definition': 'Any L-alpha-amino acid carrying an '
                                        'N-acyl substituent.',
                          'parents': ['CHEBI:58937', 'CHEBI:58937']},
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