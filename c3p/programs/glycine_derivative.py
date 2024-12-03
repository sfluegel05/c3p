"""
Classifies: CHEBI:24373 glycine derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_glycine_derivative(smiles: str):
    """
    Determines if a molecule is a glycine derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycine derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycine SMILES: NCC(=O)O
    glycine = Chem.MolFromSmiles('NCC(=O)O')
    if glycine is None:
        return False, "Error creating glycine molecule"

    # Check if the molecule contains the glycine core structure
    glycine_match = mol.HasSubstructMatch(glycine)
    if not glycine_match:
        return False, "Molecule does not contain the glycine core structure"

    # Check for modifications at the amino group or carboxy group, or replacement of any hydrogen by a heteroatom
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetDegree() > 1:
            return True, "Modification at the amino group"
        if atom.GetSymbol() == 'C' and atom.GetDegree() > 3:
            return True, "Modification at the carboxy group"
        if atom.GetSymbol() != 'C' and atom.GetSymbol() != 'H':
            return True, "Replacement of hydrogen by a heteroatom"

    return False, "No modifications found"

# Test the function with a known glycine derivative
smiles = "C(CNC(=O)CCC)(=O)O"  # butyrylglycine
print(is_glycine_derivative(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24373',
                          'name': 'glycine derivative',
                          'definition': 'A proteinogenic amino acid derivative '
                                        'resulting from reaction of glycine at '
                                        'the amino group or the carboxy group, '
                                        'or from the replacement of any '
                                        'hydrogen of glycine by a heteroatom.',
                          'parents': ['CHEBI:83811']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 7-8: malformed \\N character escape (<string>, line 1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}