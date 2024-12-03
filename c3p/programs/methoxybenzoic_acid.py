"""
Classifies: CHEBI:25238 methoxybenzoic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_methoxybenzoic_acid(smiles: str):
    """
    Determines if a molecule is a methoxybenzoic acid (benzoic acid carrying one or more methoxy substituents).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methoxybenzoic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for benzoic acid core (C1=CC=C(C=C1)C(=O)O)
    benzoic_acid = Chem.MolFromSmarts('c1ccccc1C(=O)O')
    if not mol.HasSubstructMatch(benzoic_acid):
        return False, "No benzoic acid core found"

    # Check for methoxy groups (OCH3)
    methoxy_group = Chem.MolFromSmarts('CO')
    methoxy_matches = mol.GetSubstructMatches(methoxy_group)
    if not methoxy_matches:
        return False, "No methoxy groups found"

    # Check if methoxy groups are attached to the benzene ring
    ring_info = mol.GetRingInfo()
    aromatic_rings = [ring for ring in ring_info.AtomRings() if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)]

    if not aromatic_rings:
        return False, "No aromatic rings found"

    for ring in aromatic_rings:
        if any(mol.GetAtomWithIdx(i).GetSymbol() == 'C' and mol.GetAtomWithIdx(i).GetTotalNumHs() == 0 for i in ring):
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 1:
                        return True, "Methoxybenzoic acid with methoxy substituents"
    
    return False, "Methoxy groups not attached to benzene ring"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25238',
                          'name': 'methoxybenzoic acid',
                          'definition': 'Any benzoic acid carrying one or more '
                                        'methoxy substituents.',
                          'parents': ['CHEBI:22723', 'CHEBI:51683']},
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
    'num_true_positives': 13,
    'num_false_positives': 0,
    'num_true_negatives': 16,
    'num_false_negatives': 3,
    'precision': 1.0,
    'recall': 0.8125,
    'f1': 0.896551724137931,
    'accuracy': None}