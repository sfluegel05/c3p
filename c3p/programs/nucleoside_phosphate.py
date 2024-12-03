"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate (nucleoside with one or more sugar hydroxy groups converted into mono- or poly-phosphate).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a nucleobase
    nucleobases = ['A', 'G', 'C', 'T', 'U', 'I']
    nucleobase_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in nucleobases:
            nucleobase_found = True
            break

    if not nucleobase_found:
        return False, "No nucleobase found"

    # Check for sugar moiety (ribose or deoxyribose)
    sugar_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and any(neighbor.GetSymbol() == 'C' for neighbor in atom.GetNeighbors()):
            sugar_found = True
            break

    if not sugar_found:
        return False, "No sugar moiety found"

    # Check for phosphate group(s)
    phosphate_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'P' and any(neighbor.GetSymbol() == 'O' for neighbor in atom.GetNeighbors()):
            phosphate_found = True
            break

    if not phosphate_found:
        return False, "No phosphate group found"

    return True, "Nucleoside phosphate"

# Example usage:
# smiles = "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]1O"  # Example SMILES for a nucleoside phosphate
# result, reason = is_nucleoside_phosphate(smiles)
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25608',
                          'name': 'nucleoside phosphate',
                          'definition': 'A nucleobase-containing molecular '
                                        'entity that is a nucleoside in which '
                                        'one or more of the sugar hydroxy '
                                        'groups has been converted into a '
                                        'mono- or poly-phosphate. The term '
                                        'includes both nucleotides and '
                                        'non-nucleotide nucleoside phosphates.',
                          'parents': [   'CHEBI:25703',
                                         'CHEBI:37734',
                                         'CHEBI:61120']},
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
    'num_true_positives': 163,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 0,
    'precision': 0.9878787878787879,
    'recall': 1.0,
    'f1': 0.9939024390243903,
    'accuracy': None}