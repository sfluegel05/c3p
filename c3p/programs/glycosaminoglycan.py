"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan (polysaccharide containing a substantial proportion of aminomonosaccharide residues).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify monosaccharide residues
    monosaccharide_count = 0
    amino_monosaccharide_count = 0

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            # Check for monosaccharide structure
            neighbors = atom.GetNeighbors()
            oxygen_count = sum(1 for neighbor in neighbors if neighbor.GetSymbol() == 'O')
            nitrogen_count = sum(1 for neighbor in neighbors if neighbor.GetSymbol() == 'N')
            if oxygen_count >= 2:  # Considering it as a part of sugar ring
                monosaccharide_count += 1
                if nitrogen_count >= 1:
                    amino_monosaccharide_count += 1

    if monosaccharide_count == 0:
        return False, "No monosaccharide residues found"
    
    if amino_monosaccharide_count / monosaccharide_count < 0.5:
        return False, "Aminomonosaccharide residues are not substantial"

    return True, "Contains substantial proportion of aminomonosaccharide residues"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18085',
                          'name': 'glycosaminoglycan',
                          'definition': 'Any polysaccharide containing a '
                                        'substantial proportion of '
                                        'aminomonosaccharide residues.',
                          'parents': ['CHEBI:22506']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '[19:52:55] SMILES Parse Error: syntax error while parsing: '
             'Cl/C=C\x01/CC(C(C)(C)C)CC(=O)CC=CN(C(C(CCCC1)C)=O)C\n'
             '[19:52:55] SMILES Parse Error: Failed parsing SMILES '
             "'Cl/C=C\x01/CC(C(C)(C)C)CC(=O)CC=CN(C(C(CCCC1)C)=O)C' for input: "
             "'Cl/C=C\x01/CC(C(C)(C)C)CC(=O)CC=CN(C(C(CCCC1)C)=O)C'\n"
             '[19:52:55] SMILES Parse Error: syntax error while parsing: '
             'O[C@](\\C=C\\[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])C)(C(O)(CO)C)C\n'
             '[19:52:55] SMILES Parse Error: Failed parsing SMILES '
             "'O[C@](\\C=C\\[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])C)(C(O)(CO)C)C' "
             'for input: '
             "'O[C@](\\C=C\\[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])C)(C(O)(CO)C)C'\n",
    'stdout': '',
    'num_true_positives': 3,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 26,
    'precision': 0.75,
    'recall': 0.10344827586206896,
    'f1': 0.18181818181818182,
    'accuracy': None}