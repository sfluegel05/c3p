"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole (a compound composed of two or more pyrrole units).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypyrrole, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the pyrrole substructure
    pyrrole = Chem.MolFromSmarts('c1cc[nH]c1')

    # Find all pyrrole substructures in the molecule
    pyrrole_matches = mol.GetSubstructMatches(pyrrole)

    if len(pyrrole_matches) < 2:
        return False, "Less than two pyrrole units found"

    return True, f"Found {len(pyrrole_matches)} pyrrole units"

# Example usage
smiles_examples = [
    'CCc1c(CC)c2cc3[nH]c(cc4nc(cc5[nH]c(cc1n2)c(CC)c5CC)c(CC)c4CC)c(CC)c3CC',
    'Cc1c(CCC(O)=O)c2C=c3c(\C=C\C(O)=O)c(C)c4=CC5=[N+]6C(=CC7=[N+]8C(=Cc1n2[Fe--]68n34)C(=O)[C@]7(C)CC(O)=O)C(=O)[C@]5(C)CC(O)=O',
    'C=12N3C(=CC4=[N+]5C(=CC=6N7C=8C(=C9[N+](=C(C1)[C@H]([C@@H]9CCC(O)=O)C)[Mg-2]753)CC(C8C6CC)=O)C(=C4C)CCC)C(=C2C)C=C'
]

for smiles in smiles_examples:
    print(is_polypyrrole(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38077',
                          'name': 'polypyrrole',
                          'definition': 'A compound composed of two or more '
                                        'pyrrole units.',
                          'parents': ['CHEBI:38101']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[00:11:10] SMILES Parse Error: syntax error while parsing: '
             'C[C@]1(CC(O)=O)[C@H](CCC(O)=O)C2=C/c3[nH]c(Cc4[nH]c(c(CC(O)=O)c4CCC(O)=O)[C@](C)(O)[C@@]45N/C(=C\\C1=N\x02)[C@@H](CCC(O)=O)[C@]4(C)CC(=O)O5)c(CCC(O)=O)c3CC(O)=O\n'
             '[00:11:10] SMILES Parse Error: Failed parsing SMILES '
             "'C[C@]1(CC(O)=O)[C@H](CCC(O)=O)C2=C/c3[nH]c(Cc4[nH]c(c(CC(O)=O)c4CCC(O)=O)[C@](C)(O)[C@@]45N/C(=C\\C1=N\x02)[C@@H](CCC(O)=O)[C@]4(C)CC(=O)O5)c(CCC(O)=O)c3CC(O)=O' "
             'for input: '
             "'C[C@]1(CC(O)=O)[C@H](CCC(O)=O)C2=C/c3[nH]c(Cc4[nH]c(c(CC(O)=O)c4CCC(O)=O)[C@](C)(O)[C@@]45N/C(=C\\C1=N\x02)[C@@H](CCC(O)=O)[C@]4(C)CC(=O)O5)c(CCC(O)=O)c3CC(O)=O'\n",
    'stdout': "(True, 'Found 2 pyrrole units')\n"
              "(False, 'Less than two pyrrole units found')\n"
              "(False, 'Less than two pyrrole units found')\n",
    'num_true_positives': 14,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 23,
    'precision': 1.0,
    'recall': 0.3783783783783784,
    'f1': 0.5490196078431372,
    'accuracy': None}