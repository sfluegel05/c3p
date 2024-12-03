"""
Classifies: CHEBI:22221 acyl group
"""
from rdkit import Chem

def is_acyl_group(smiles: str):
    """
    Determines if a molecule is an acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl group, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carbonyl group (C=O)
    carbonyl = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() == 2.0:
                    carbonyl = True
                    break
        if carbonyl:
            break

    if not carbonyl:
        return False, "No carbonyl group (C=O) found"

    # Check for the presence of a carbon atom with a single bond to a hydrogen or heteroatom (excluding oxygen in carbonyl)
    acyl_bond = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            if atom.GetTotalNumHs() == 1 or any(neighbor.GetAtomicNum() not in (6, 8) for neighbor in atom.GetNeighbors()):
                acyl_bond = True
                break

    if not acyl_bond:
        return False, "No acyl bond found (carbon with single bond to hydrogen or heteroatom)"

    return True, "Valid acyl group"

# Example usage
print(is_acyl_group("CCCCCCCCCC\C=C/CCC(-*)=O"))  # (4Z)-hexadecenoyl group
print(is_acyl_group("C1=CC(=CC=C1)C(*)=O"))  # benzoyl group
print(is_acyl_group("C(CNC(CCCCCCCCCCCCCCC)=O)(=O)*"))  # N-hexadecanoylglycyl group
print(is_acyl_group("CC"))  # Invalid example, should return False, "No carbonyl group (C=O) found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22221',
                          'name': 'acyl group',
                          'definition': 'An organic group formed by removing '
                                        'one or more hydroxy groups from an '
                                        'oxoacid that has the general '
                                        'structure RkE(=O)l(OH)m (l =/= 0). '
                                        'Although the term is almost always '
                                        'applied to organic compounds, with '
                                        'carboxylic acid as the oxoacid, acyl '
                                        'groups can in principle be derived '
                                        'from other types of acids such as '
                                        'sulfonic acids or phosphonic acids.',
                          'parents': ['CHEBI:33247']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Valid acyl group')\n"
              "(True, 'Valid acyl group')\n"
              "(True, 'Valid acyl group')\n"
              "(False, 'No carbonyl group (C=O) found')\n",
    'num_true_positives': 33,
    'num_false_positives': 19,
    'num_true_negatives': 1,
    'num_false_negatives': 1,
    'precision': 0.6346153846153846,
    'recall': 0.9705882352941176,
    'f1': 0.7674418604651162,
    'accuracy': None}