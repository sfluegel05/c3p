"""
Classifies: CHEBI:35317 glycosylgalactose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_glycosylgalactose(smiles: str):
    """
    Determines if a molecule is a glycosylgalactose (a disaccharide which has a galactose residue at the reducing end).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosylgalactose, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMILES of galactose
    galactose_smiles = 'C(C1C(C(C(C(O1)O)O)O)O)O'
    galactose_mol = Chem.MolFromSmiles(galactose_smiles)
    if galactose_mol is None:
        return False, "Error in defining galactose structure"

    # Check if the molecule contains galactose as a substructure
    if not mol.HasSubstructMatch(galactose_mol):
        return False, "No galactose residue found"

    # Check if the galactose residue is at the reducing end
    reducing_end_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetDegree() == 1:
            neighbor = atom.GetNeighbors()[0]
            if neighbor.HasSubstructMatch(galactose_mol):
                reducing_end_found = True
                break

    if not reducing_end_found:
        return False, "Galactose residue is not at the reducing end"

    # Check if the molecule is a disaccharide
    num_sugars = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetDegree() == 2:
            num_sugars += 1

    if num_sugars != 2:
        return False, "Molecule is not a disaccharide"

    return True, "Molecule is a glycosylgalactose"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35317',
                          'name': 'glycosylgalactose',
                          'definition': 'A disaccharide which has a galactose '
                                        'residue at the reducing end.',
                          'parents': ['CHEBI:36233']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "'Atom' object has no attribute 'HasSubstructMatch'",
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}