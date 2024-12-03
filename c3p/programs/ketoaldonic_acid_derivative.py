"""
Classifies: CHEBI:63394 ketoaldonic acid derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_ketoaldonic_acid_derivative(smiles: str):
    """
    Determines if a molecule is a ketoaldonic acid derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ketoaldonic acid derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a keto group (C=O)
    has_keto_group = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            for bond in atom.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetOtherAtom(atom).GetAtomicNum() == 8:
                    has_keto_group = True
                    break
        if has_keto_group:
            break

    if not has_keto_group:
        return False, "No keto group (C=O) found"

    # Check for the presence of a carboxylic acid group (C(=O)O)
    has_carboxylic_acid_group = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            neighbors = [neighbor.GetAtomicNum() for neighbor in atom.GetNeighbors()]
            if neighbors.count(8) == 2:  # Two oxygens (one double-bonded, one single-bonded)
                for bond in atom.GetBonds():
                    if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetOtherAtom(atom).GetAtomicNum() == 8:
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetAtomicNum() == 8 and neighbor.GetBonds()[0].GetBondType() == Chem.BondType.SINGLE:
                                has_carboxylic_acid_group = True
                                break
        if has_carboxylic_acid_group:
            break

    if not has_carboxylic_acid_group:
        return False, "No carboxylic acid group (C(=O)O) found"

    # Check for the presence of a monosaccharide backbone
    ring_info = mol.GetRingInfo()
    if not ring_info.AtomRings():
        return False, "No ring structure found"
    
    # Check for a ring with oxygen atoms, indicating a sugar-like structure
    has_oxygen_in_ring = False
    for ring in ring_info.AtomRings():
        if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring):
            has_oxygen_in_ring = True
            break

    if not has_oxygen_in_ring:
        return False, "No oxygen-containing ring found, indicating no monosaccharide backbone"

    return True, "Ketoaldonic acid derivative found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:63394',
                          'name': 'ketoaldonic acid derivative',
                          'definition': 'A monosaccharide derivative that is '
                                        'formally obtained from a ketoaldonic '
                                        'acid.',
                          'parents': ['CHEBI:63367']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 119,
    'num_false_positives': 5,
    'num_true_negatives': 15,
    'num_false_negatives': 1,
    'precision': 0.9596774193548387,
    'recall': 0.9916666666666667,
    'f1': 0.9754098360655739,
    'accuracy': None}