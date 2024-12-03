"""
Classifies: CHEBI:23643 depsipeptide
"""
from rdkit import Chem

def is_depsipeptide(smiles: str):
    """
    Determines if a molecule is a depsipeptide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a depsipeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of both amide and ester bonds
    has_amide = False
    has_ester = False

    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 8:  # C=O
                if atom1.GetTotalDegree() == 3 and atom2.GetTotalDegree() == 1:
                    for neighbor in atom1.GetNeighbors():
                        if neighbor.GetAtomicNum() == 7:  # N
                            has_amide = True
                        elif neighbor.GetAtomicNum() == 8:  # O
                            has_ester = True
                elif atom2.GetTotalDegree() == 3 and atom1.GetTotalDegree() == 1:
                    for neighbor in atom2.GetNeighbors():
                        if neighbor.GetAtomicNum() == 7:  # N
                            has_amide = True
                        elif neighbor.GetAtomicNum() == 8:  # O
                            has_ester = True

    if has_amide and has_ester:
        return True, "Molecule is a depsipeptide"
    else:
        return False, "Molecule does not contain both amide and ester bonds"

# Example usage:
# smiles = "O=C1OCC(O)CC(NC(=O)C(NC(=O)C(NC(=O)[C@H]2NCCC2)CCC(=O)N)CC3=CC=C(O)C=C3)C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(NC(C(NC(C(NC1C(CC)C)=O)CC4=CC=C(O)C=C4)=O)C)=O)CCCN)=O)CC5=CC=C(O)C=C5)C(O)C)CCC(=O)O)C(C)C"
# result, reason = is_depsipeptide(smiles)
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23643',
                          'name': 'depsipeptide',
                          'definition': 'A natural or synthetic compound '
                                        'having a sequence of amino and '
                                        'hydroxy carboxylic acid residues '
                                        '(usually alpha-amino and '
                                        'alpha-hydroxy acids), commonly but '
                                        'not necessarily regularly '
                                        'alternating.',
                          'parents': ['CHEBI:16670']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 67-68: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}