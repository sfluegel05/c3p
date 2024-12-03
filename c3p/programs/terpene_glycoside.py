"""
Classifies: CHEBI:61777 terpene glycoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_terpene_glycoside(smiles: str):
    """
    Determines if a molecule is a terpene glycoside.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a terpene glycoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycosidic bonds
    glycosidic_bond = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and len(atom.GetNeighbors()) == 2:
            neighbors = atom.GetNeighbors()
            if any(neigh.GetSymbol() == 'C' and neigh.GetDegree() > 1 for neigh in neighbors):
                glycosidic_bond = True
                break

    if not glycosidic_bond:
        return False, "No glycosidic bond found"

    # Check for terpenoid backbone
    terpenoid_backbone = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            # Check for isoprene units (C5H8)
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 4:
                c_count = sum(1 for neigh in neighbors if neigh.GetSymbol() == 'C')
                if c_count == 2:
                    terpenoid_backbone = True
                    break

    if not terpenoid_backbone:
        return False, "No terpenoid backbone found"

    return True, "Molecule is a terpene glycoside"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61777',
                          'name': 'terpene glycoside',
                          'definition': 'A terpenoid in which one or more '
                                        'hydroxy functions are glycosylated.',
                          'parents': ['CHEBI:24400', 'CHEBI:26873']},
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
    'num_true_positives': 11,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 87,
    'precision': 1.0,
    'recall': 0.11224489795918367,
    'f1': 0.20183486238532108,
    'accuracy': None}