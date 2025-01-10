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
    if not amino_acid_matches:
        return False, "No L-alpha-amino acid backbone found"
    
    # Get the nitrogen atom of the amino group
    amino_nitrogen = amino_acid_matches[0][2]

    # Check if the nitrogen has exactly one acyl group attached
    nitrogen_atom = mol.GetAtomWithIdx(amino_nitrogen)
    if nitrogen_atom.GetDegree() != 2:  # Should be connected to H and acyl group
        return False, "Amino nitrogen has incorrect number of substituents"

    # Check for N-acyl substituent: must be directly attached to the amino nitrogen
    acyl_group = None
    for neighbor in nitrogen_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # Carbon
            for bond in neighbor.GetBonds():
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and \
                   bond.GetBeginAtom().GetAtomicNum() == 8 and \
                   bond.GetEndAtom().GetAtomicNum() == 6:
                    acyl_group = neighbor
                    break
    
    if not acyl_group:
        return False, "No N-acyl substituent found attached to amino nitrogen"

    # Verify the acyl group is not part of a peptide bond
    # Should only have one connection to the nitrogen and one to another carbon
    if acyl_group.GetDegree() != 2:
        return False, "Acyl group is part of a peptide bond"

    # Verify stereochemistry at alpha carbon
    alpha_carbon = amino_acid_matches[0][0]
    if not mol.GetAtomWithIdx(alpha_carbon).GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
        return False, "Incorrect stereochemistry at alpha carbon"

    return True, "Contains L-alpha-amino acid backbone with N-acyl substituent"


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