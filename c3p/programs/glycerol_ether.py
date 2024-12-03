"""
Classifies: CHEBI:24353 glycerol ether
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glycerol_ether(smiles: str):
    """
    Determines if a molecule is a glycerol ether.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerol ether, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMILES pattern for glyceryl ether
    glycerol_ether_pattern = Chem.MolFromSmarts("OCC(O)CO")

    # Check if the molecule contains the glyceryl ether pattern
    if not mol.HasSubstructMatch(glycerol_ether_pattern):
        return False, "No glyceryl ether pattern found"

    # Check if there is at least one ether linkage with glyceryl as O-substituent
    glyceryl_ether_matches = mol.GetSubstructMatches(glycerol_ether_pattern)
    for match in glyceryl_ether_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'O':
                neighbors = atom.GetNeighbors()
                for neighbor in neighbors:
                    if neighbor.GetSymbol() == 'C':
                        if any(n.GetSymbol() == 'O' for n in neighbor.GetNeighbors()):
                            return True, "Glycerol ether found"

    return False, "No ether linkage with glyceryl as O-substituent found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24353',
                          'name': 'glycerol ether',
                          'definition': 'Any ether having glyceryl as at least '
                                        'one of the O-substituents.',
                          'parents': ['CHEBI:25698']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[20:28:19] SMILES Parse Error: syntax error while parsing: '
             'O[C@@H]1C/C(=C\\C=C\x02/[C@]3([C@@]([C@](CC3)([C@H](C)C=C[C@@H](C(C)C)C)[H])(CCC2)C)[H])/C([C@@H](O)C1)=C\n'
             '[20:28:19] SMILES Parse Error: Failed parsing SMILES '
             "'O[C@@H]1C/C(=C\\C=C\x02/[C@]3([C@@]([C@](CC3)([C@H](C)C=C[C@@H](C(C)C)C)[H])(CCC2)C)[H])/C([C@@H](O)C1)=C' "
             'for input: '
             "'O[C@@H]1C/C(=C\\C=C\x02/[C@]3([C@@]([C@](CC3)([C@H](C)C=C[C@@H](C(C)C)C)[H])(CCC2)C)[H])/C([C@@H](O)C1)=C'\n",
    'stdout': '',
    'num_true_positives': 54,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 4,
    'precision': 0.9642857142857143,
    'recall': 0.9310344827586207,
    'f1': 0.9473684210526316,
    'accuracy': None}