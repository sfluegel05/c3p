"""
Classifies: CHEBI:46668 amino acid ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an amino acid ester (carboxylic ester derivative of an amino acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino acid ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ester functional group (COOC)
    ester_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 3:
                oxygens = [n for n in neighbors if n.GetSymbol() == 'O']
                if len(oxygens) == 2:
                    ester_found = True
                    break

    if not ester_found:
        return False, "No ester functional group found"

    # Check for amino acid structure
    amino_acid_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            neighbors = atom.GetNeighbors()
            carbons = [n for n in neighbors if n.GetSymbol() == 'C']
            if len(carbons) >= 1:
                for carbon in carbons:
                    carboxyl_group = any(
                        neighbor.GetSymbol() == 'O' and neighbor.GetTotalDegree() == 1
                        for neighbor in carbon.GetNeighbors()
                    )
                    if carboxyl_group:
                        amino_acid_found = True
                        break
            if amino_acid_found:
                break

    if not amino_acid_found:
        return False, "No amino acid structure found"

    return True, "Molecule is an amino acid ester"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:46668',
                          'name': 'amino acid ester',
                          'definition': 'Any  carboxylic ester derivative of '
                                        'an amino acid.',
                          'parents': ['CHEBI:33308', 'CHEBI:83821']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 16,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 14,
    'precision': 1.0,
    'recall': 0.5333333333333333,
    'f1': 0.6956521739130436,
    'accuracy': None}