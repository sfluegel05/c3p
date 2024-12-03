"""
Classifies: CHEBI:52395 oxyketone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oxyketone(smiles: str):
    """
    Determines if a molecule is an oxyketone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxyketone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a ketone group (C=O)
    ketone = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:  # Oxygen
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        ketone = True
                        break
            if ketone:
                break

    if not ketone:
        return False, "No ketone group (C=O) found"

    # Check for the presence of an oxy group (-O-) in the R groups
    oxy_group = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 2 and all(neighbor.GetAtomicNum() == 6 for neighbor in neighbors):
                oxy_group = True
                break

    if not oxy_group:
        return False, "No oxy group (-O-) found in the R groups"

    return True, "Oxyketone found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:52395',
                          'name': 'oxyketone',
                          'definition': 'A compound with the general formula '
                                        'R2C=O (R=/=H) where one or more of '
                                        'the R groups contains an oxy (-O-) '
                                        'group.',
                          'parents': ['CHEBI:17087']},
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
    'num_true_positives': 21,
    'num_false_positives': 11,
    'num_true_negatives': 9,
    'num_false_negatives': 30,
    'precision': 0.65625,
    'recall': 0.4117647058823529,
    'f1': 0.5060240963855422,
    'accuracy': None}