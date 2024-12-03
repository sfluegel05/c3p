"""
Classifies: CHEBI:24399 glycosaminyl group
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_glycosaminyl_group(smiles: str):
    """
    Determines if a molecule is a glycosaminyl group.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a glycosaminyl group, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycosamine structure
    glycosamine_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    for nneighbor in neighbor.GetNeighbors():
                        if nneighbor.GetSymbol() == 'O':
                            glycosamine_found = True
                            break
                    if glycosamine_found:
                        break
            if glycosamine_found:
                break

    if not glycosamine_found:
        return False, "No glycosamine structure found"

    # Check for glycosyl group
    glycosyl_group_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    for nneighbor in neighbor.GetNeighbors():
                        if nneighbor.GetSymbol() == 'O':
                            glycosyl_group_found = True
                            break
                    if glycosyl_group_found:
                        break
            if glycosyl_group_found:
                break

    if not glycosyl_group_found:
        return False, "No glycosyl group found"

    # Check for removal of hydroxy group from the hemiacetal function
    hemiacetal_function_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            oxy_count = 0
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O':
                    oxy_count += 1
            if oxy_count == 1:
                hemiacetal_function_found = True
                break

    if not hemiacetal_function_found:
        return False, "No hemiacetal function found"

    return True, "Molecule is a glycosaminyl group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24399',
                          'name': 'glycosaminyl group',
                          'definition': 'A glycosyl group obtained by removing '
                                        'the hydroxy group from the hemiacetal '
                                        'function of a glycosamine and, by '
                                        'extension, of a lower oligosaccharide '
                                        'having a glycosamine at the reducing '
                                        'end.',
                          'parents': ['CHEBI:24403']},
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
    'num_false_positives': 4,
    'num_true_negatives': 7,
    'num_false_negatives': 0,
    'precision': 0.7333333333333333,
    'recall': 1.0,
    'f1': 0.846153846153846,
    'accuracy': None}