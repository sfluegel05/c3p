"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine (an alkylamine in which the alkyl group is substituted by an aromatic group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule contains an amine group (NH2, NH, or N)
    amine_groups = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetDegree() > 0]
    if not amine_groups:
        return False, "No amine group found"

    # Check if the amine group is attached to an alkyl chain
    for amine in amine_groups:
        neighbors = amine.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:  # Carbon
                # Check if the carbon is part of an alkyl chain
                if any(alkyl_neighbor.GetAtomicNum() == 6 for alkyl_neighbor in neighbor.GetNeighbors() if alkyl_neighbor.GetIdx() != amine.GetIdx()):
                    # Check if the alkyl chain is substituted by an aromatic group
                    for alkyl_neighbor in neighbor.GetNeighbors():
                        if alkyl_neighbor.GetAtomicNum() == 6 and alkyl_neighbor.GetIdx() != amine.GetIdx():
                            if any(aromatic_neighbor.GetIsAromatic() for aromatic_neighbor in alkyl_neighbor.GetNeighbors()):
                                return True, "Aralkylamine found"
    
    return False, "No aralkylamine structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18000',
                          'name': 'aralkylamine',
                          'definition': 'An alkylamine in which the alkyl '
                                        'group is substituted by an aromatic '
                                        'group.',
                          'parents': ['CHEBI:22331', 'CHEBI:64365']},
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
    'num_true_positives': 12,
    'num_false_positives': 13,
    'num_true_negatives': 0,
    'num_false_negatives': 1,
    'precision': 0.48,
    'recall': 0.9230769230769231,
    'f1': 0.631578947368421,
    'accuracy': None}