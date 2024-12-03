"""
Classifies: CHEBI:26273 proline derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_proline_derivative(smiles: str):
    """
    Determines if a molecule is a proline derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proline derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMILES string for proline
    proline_smiles = "C1CC(NC1)C(=O)O"
    proline = Chem.MolFromSmiles(proline_smiles)
    
    if proline is None:
        return False, "Error in defining proline structure"

    # Check if the molecule contains the proline substructure
    if not mol.HasSubstructMatch(proline):
        return False, "Molecule does not contain proline substructure"

    # Check for modifications at the amino group or carboxy group
    proline_amino_group = Chem.MolFromSmarts("C1CC(NC1)C(=O)O")
    proline_carboxy_group = Chem.MolFromSmarts("C1CC(NC1)C(=O)O")
    
    if not (mol.HasSubstructMatch(proline_amino_group) or mol.HasSubstructMatch(proline_carboxy_group)):
        return False, "Molecule does not have modifications at the amino or carboxy group"

    # Check for replacement of any hydrogen by a heteroatom
    heteroatoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6)]
    if not heteroatoms:
        return False, "No heteroatoms found in the molecule"

    # Exclude peptides containing proline residues
    peptide_bond = Chem.MolFromSmarts("N[C@@H](C(=O))")
    if mol.HasSubstructMatch(peptide_bond):
        return False, "Molecule contains peptide bonds"

    return True, "Molecule is a proline derivative"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26273',
                          'name': 'proline derivative',
                          'definition': 'An amino acid derivative resulting '
                                        'from reaction of proline at the amino '
                                        'group or the carboxy group, or from '
                                        'the replacement of any hydrogen of '
                                        'proline by a heteroatom. The '
                                        'definition normally excludes peptides '
                                        'containing proline residues.',
                          'parents': ['CHEBI:83821']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 15,
    'num_false_negatives': 15,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}