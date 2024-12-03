"""
Classifies: CHEBI:33710 alpha-amino-acid residue
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_alpha_amino_acid_residue(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid residue.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid residue, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of an alpha-amino acid backbone
    alpha_amino_acid_backbone = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetTotalDegree() >= 1:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and neighbor.GetTotalDegree() == 4:
                    for n_neighbor in neighbor.GetNeighbors():
                        if n_neighbor.GetSymbol() == 'C' and n_neighbor.GetTotalDegree() >= 2 and n_neighbor.GetIdx() != atom.GetIdx():
                            alpha_amino_acid_backbone = True
                            break
                    if alpha_amino_acid_backbone:
                        break
            if alpha_amino_acid_backbone:
                break

    if not alpha_amino_acid_backbone:
        return False, "No alpha-amino acid backbone found"

    # Check for the presence of a carbonyl group attached to the alpha carbon
    carbonyl_group = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() >= 2:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O' and (neighbor.GetTotalDegree() == 1 or neighbor.GetTotalDegree() == 2):
                    carbonyl_group = True
                    break
            if carbonyl_group:
                break

    if not carbonyl_group:
        return False, "No carbonyl group attached to the alpha carbon found"

    return True, "Alpha-amino acid residue found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33710',
                          'name': 'alpha-amino-acid residue',
                          'definition': 'An amino-acid residue derived from an '
                                        'alpha-amino acid.',
                          'parents': ['CHEBI:33708']},
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
    'num_true_positives': 52,
    'num_false_positives': 10,
    'num_true_negatives': 10,
    'num_false_negatives': 0,
    'precision': 0.8387096774193549,
    'recall': 1.0,
    'f1': 0.9122807017543859,
    'accuracy': None}