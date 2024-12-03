"""
Classifies: CHEBI:38834 benzofurochromene
"""
from rdkit import Chem

def is_benzofurochromene(smiles: str):
    """
    Determines if a molecule is a benzofurochromene or its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a benzofurochromene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for benzofurochromene core
    benzofurochromene_smarts = 'O1C2=C(C=3C1=CC=4OCOC4C3)COC5=C2C=CC=6C=CC=CC6=C5'

    # Create a molecule object from the SMARTS pattern
    benzofurochromene_mol = Chem.MolFromSmarts(benzofurochromene_smarts)
    if benzofurochromene_mol is None:
        return False, "Invalid SMARTS pattern for benzofurochromene"

    # Check for the presence of the benzofurochromene core in the molecule
    if mol.HasSubstructMatch(benzofurochromene_mol):
        return True, "Molecule is a benzofurochromene or its substituted derivative"
    else:
        return False, "Molecule does not contain the benzofurochromene core"

# Example usage
smiles_examples = [
    "O1C2=C(C=3C1=CC(O)=C(O)C3)C(OC=4C2=C(O)C=C(O)C4)=O",
    "O1C2C(O)(C3=C1C=C(OC)C=C3)COC4=C2C=C5C(OC(=C5)C(O)(C)C)=C4",
    "O1C2C(C3=C1C=C(OC)C=C3)COC=4C2=CC(OC)=C(O)C4",
    "[H][C@@]12Oc3cc(O)ccc3[C@]1(O)COc1cc(O)c(CC=C(C)C)cc21",
    "O1[C@@]2([C@](C3=C1C=C(O)C=C3)(COC4=C2C=CC(OC)=C4)[H])[H]",
    "COc1c(CC=C(C)C)c(O)cc2OCc3c(oc4cc(O)c(CC=C(C)C)cc34)-c12",
    "COc1ccc2[C@@H]3COc4cc(O)ccc4[C@@H]3Oc2c1OC",
    "O1C2C(C3=C1C(OC)=C(OC)C=C3)COC4=C2C=CC(OC)=C4O",
    "COc1cc(O)c2c(c1)oc(=O)c1c3cc(O)c(O)cc3oc21",
    "O1C2C(C3=C1C(C[C@H](O)C(C)=C)=C(O)C=C3)COC4=C2C=CC(O)=C4",
    "COc1ccc2C3COc4cc(O)ccc4C3Oc2c1",
    "O1C2=C(C=3C1=CC=4OCOC4C3)COC5=C2C=CC(O)=C5",
    "O1C(C(O)C(O)C2=C1C=C3OC(=O)C4=C(OC5=C4C=CC(O)=C5)C3=C2)(C)C",
    "[H][C@@]12COc3cc(OC)ccc3[C@]1([H])Oc1c(CC=C(C)C)c(OC)c(O)cc21",
    "O1C2C(C3=C1C=C(O)C=C3)COC4=C2C=CC(OC)=C4O",
]

for smiles in smiles_examples:
    result, reason = is_benzofurochromene(smiles)
    print(f"SMILES: {smiles}, Result: {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38834',
                          'name': 'benzofurochromene',
                          'definition': 'An organic heteropolycyclic compound '
                                        'whose skeleton consists of a chromene '
                                        'ring fused onto a 1-benzofuran ring, '
                                        'together with their substituted '
                                        'derivatives.',
                          'parents': ['CHEBI:38166']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 21-22: malformed \\N character escape (<string>, line '
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