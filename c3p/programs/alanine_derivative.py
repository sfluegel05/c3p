"""
Classifies: CHEBI:22278 alanine derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_alanine_derivative(smiles: str):
    """
    Determines if a molecule is an alanine derivative, defined as an amino acid derivative resulting from reaction of alanine at the amino group or the carboxy group, or from the replacement of any hydrogen of alanine by a heteroatom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alanine derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the alanine structure
    alanine_smiles = "CC(C(=O)O)N"
    alanine_mol = Chem.MolFromSmiles(alanine_smiles)

    # Check for substructure match
    if not mol.HasSubstructMatch(alanine_mol):
        return False, "The molecule does not contain the alanine core structure"

    # Check for reactions at the amino group or carboxy group, or replacement of any hydrogen by a heteroatom
    alanine_atoms = [atom.GetIdx() for atom in alanine_mol.GetAtoms()]
    mol_atoms = [mol.GetAtomWithIdx(idx) for idx in range(mol.GetNumAtoms())]
    
    for idx in alanine_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() != 'C' and atom.GetSymbol() != 'N' and atom.GetSymbol() != 'O':
            return True, "The molecule is an alanine derivative with heteroatom substitution"
    
    for bond in alanine_mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if (mol_atoms[begin_idx].GetSymbol() != alanine_mol.GetAtomWithIdx(begin_idx).GetSymbol() or
            mol_atoms[end_idx].GetSymbol() != alanine_mol.GetAtomWithIdx(end_idx).GetSymbol()):
            return True, "The molecule is an alanine derivative with reaction at the amino or carboxy group"
    
    return False, "The molecule does not meet the criteria for being an alanine derivative"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22278',
                          'name': 'alanine derivative',
                          'definition': 'An amino acid derivative resulting '
                                        'from reaction of alanine at the amino '
                                        'group or the carboxy group, or from '
                                        'the replacement of any hydrogen of '
                                        'alanine by a heteroatom. The '
                                        'definition normally excludes peptides '
                                        'containing alanine residues.',
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
    'num_true_positives': 11,
    'num_false_positives': 1,
    'num_true_negatives': 12,
    'num_false_negatives': 2,
    'precision': 0.9166666666666666,
    'recall': 0.8461538461538461,
    'f1': 0.8799999999999999,
    'accuracy': None}