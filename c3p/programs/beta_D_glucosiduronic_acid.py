"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define beta-D-glucuronic acid fragment
    beta_D_glucuronic_acid_smiles = "O[C@@H]1[C@@H](O)[C@H](O)[C@@H](O[C@H]1C(O)=O)O"
    beta_D_glucuronic_acid_mol = Chem.MolFromSmiles(beta_D_glucuronic_acid_smiles)
    if beta_D_glucuronic_acid_mol is None:
        return False, "Invalid beta-D-glucuronic acid SMILES string"

    # Check if the molecule contains the beta-D-glucuronic acid fragment
    if not mol.HasSubstructMatch(beta_D_glucuronic_acid_mol):
        return False, "Molecule does not contain beta-D-glucuronic acid fragment"
    
    # Check for glycosidic bond
    matches = mol.GetSubstructMatches(beta_D_glucuronic_acid_mol)
    for match in matches:
        beta_D_glucuronic_acid_atoms = set(match)
        for atom_idx in beta_D_glucuronic_acid_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in beta_D_glucuronic_acid_atoms:
                    return True, "Molecule is a beta-D-glucosiduronic acid"
    
    return False, "No glycosidic bond found"

# Example usage
smiles_examples = [
    "O1[C@@H]([C@H]([C@@H]([C@H]([C@@H]1OC=2C3=C(C=4C=CC=CC4C2)C=CC=C3)O)O)O)C(O)=O",
    "C1[C@]2([C@]3([C@@](C4=C(C=C(O)C(=C4)O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)C(=O)O)O)O)O)CC3)(CC[C@@]2(C(=O)C1)C)[H])[H])[H]",
    "O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OC2=C(OC)C=C(C=C2)C=3OC4=C(C(=O)C3OC)C(O)=C(OC)C(O)=C4)C(O)=O",
    "C=1C(=CC(=C2C1OC(C(C2=O)O)C3=CC(=C(C=C3)O)OC4C(C(C(C(O4)C(O)=O)O)O)O)O)O",
    "O1C(C(O)C(O)C(O)C1OC2=C(OC)C=C(C=C2)C=3OC=4C(C(=O)C3OC)=C(O)C=5OCOC5C4)C(O)=O",
    "CC\C(c1ccccc1)=C(/c1ccccc1)c1ccc(OCC[N+](C)(C)[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)C(O)=O)cc1",
    "C=1(C=C(C=CC1)CC2CCC(O2)=O)O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)C(O)=O)O)O)O",
    "O1C(OC=2C=C(C3OC=4C(CC3O)=C(O)C=C(O)C4)C=CC2OC)C(O)C(O)C(O)C1C(O)=O",
    "O[C@@H]1[C@@H](O)[C@H](Oc2cc(O)cc(CCc3ccc(O)cc3)c2)O[C@@H]([C@H]1O)C(O)=O",
    "O[C@@H]1[C@@H](O)[C@@H](O[C@@H]([C@H]1O)C(O)=O)OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
    "CNC1=NC2=C(C(C)=C(O[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3O)C(O)=O)C(C)=C2S1)CC4=CN=CC=C4",
    "O[C@@H]1[C@@H](O)[C@@H](O[C@@H]([C@H]1O)C(O)=O)Oc1cc(O)c2C(=O)C[C@@H](Oc2c1)c1ccc(O)c(O)c1",
    "[C@@H]1([C@@H]([C@@H](O)[C@@H]([C@H](O1)C(O)=O)O)O)OC=2C(=CC=3C=C(C=CC3C2)Br)C(=O)NC=4C(=CC=CC4)OC"
]

for smiles in smiles_examples:
    result, reason = is_beta_D_glucosiduronic_acid(smiles)
    print(f"SMILES: {smiles}\nResult: {result}, Reason: {reason}\n")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:15341',
                          'name': 'beta-D-glucosiduronic acid',
                          'definition': 'A glucosiduronic acid resulting from '
                                        'the formal condensation of any '
                                        'substance with beta-D-glucuronic acid '
                                        'to form a glycosidic bond.',
                          'parents': ['CHEBI:24302']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: '
              'O1[C@@H]([C@H]([C@@H]([C@H]([C@@H]1OC=2C3=C(C=4C=CC=CC4C2)C=CC=C3)O)O)O)C(O)=O\n'
              'Result: True, Reason: Molecule is a beta-D-glucosiduronic acid\n'
              '\n'
              'SMILES: '
              'C1[C@]2([C@]3([C@@](C4=C(C=C(O)C(=C4)O[C@H]5[C@@H]([C@H]([C@@H]([C@H](O5)C(=O)O)O)O)O)CC3)(CC[C@@]2(C(=O)C1)C)[H])[H])[H]\n'
              'Result: True, Reason: Molecule is a beta-D-glucosiduronic acid\n'
              '\n'
              'SMILES: '
              'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OC2=C(OC)C=C(C=C2)C=3OC4=C(C(=O)C3OC)C(O)=C(OC)C(O)=C4)C(O)=O\n'
              'Result: True, Reason: Molecule is a beta-D-glucosiduronic acid\n'
              '\n'
              'SMILES: '
              'C=1C(=CC(=C2C1OC(C(C2=O)O)C3=CC(=C(C=C3)O)OC4C(C(C(C(O4)C(O)=O)O)O)O)O)O\n'
              'Result: True, Reason: Molecule is a beta-D-glucosiduronic acid\n'
              '\n'
              'SMILES: '
              'O1C(C(O)C(O)C(O)C1OC2=C(OC)C=C(C=C2)C=3OC=4C(C(=O)C3OC)=C(O)C=5OCOC5C4)C(O)=O\n'
              'Result: True, Reason: Molecule is a beta-D-glucosiduronic acid\n'
              '\n'
              'SMILES: '
              'CC\\C(c1ccccc1)=C(/c1ccccc1)c1ccc(OCC[N+](C)(C)[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)C(O)=O)cc1\n'
              'Result: False, Reason: Molecule does not contain '
              'beta-D-glucuronic acid fragment\n'
              '\n'
              'SMILES: '
              'C=1(C=C(C=CC1)CC2CCC(O2)=O)O[C@H]3[C@@H]([C@H]([C@@H]([C@H](O3)C(O)=O)O)O)O\n'
              'Result: True, Reason: Molecule is a beta-D-glucosiduronic acid\n'
              '\n'
              'SMILES: '
              'O1C(OC=2C=C(C3OC=4C(CC3O)=C(O)C=C(O)C4)C=CC2OC)C(O)C(O)C(O)C1C(O)=O\n'
              'Result: True, Reason: Molecule is a beta-D-glucosiduronic acid\n'
              '\n'
              'SMILES: '
              'O[C@@H]1[C@@H](O)[C@H](Oc2cc(O)cc(CCc3ccc(O)cc3)c2)O[C@@H]([C@H]1O)C(O)=O\n'
              'Result: True, Reason: Molecule is a beta-D-glucosiduronic acid\n'
              '\n'
              'SMILES: '
              'O[C@@H]1[C@@H](O)[C@@H](O[C@@H]([C@H]1O)C(O)=O)OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl\n'
              'Result: True, Reason: Molecule is a beta-D-glucosiduronic acid\n'
              '\n'
              'SMILES: '
              'CNC1=NC2=C(C(C)=C(O[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3O)C(O)=O)C(C)=C2S1)CC4=CN=CC=C4\n'
              'Result: True, Reason: Molecule is a beta-D-glucosiduronic acid\n'
              '\n'
              'SMILES: '
              'O[C@@H]1[C@@H](O)[C@@H](O[C@@H]([C@H]1O)C(O)=O)Oc1cc(O)c2C(=O)C[C@@H](Oc2c1)c1ccc(O)c(O)c1\n'
              'Result: True, Reason: Molecule is a beta-D-glucosiduronic acid\n'
              '\n'
              'SMILES: '
              '[C@@H]1([C@@H]([C@@H](O)[C@@H]([C@H](O1)C(O)=O)O)O)OC=2C(=CC=3C=C(C=CC3C2)Br)C(=O)NC=4C(=CC=CC4)OC\n'
              'Result: True, Reason: Molecule is a beta-D-glucosiduronic acid\n'
              '\n',
    'num_true_positives': 12,
    'num_false_positives': 13,
    'num_true_negatives': 0,
    'num_false_negatives': 1,
    'precision': 0.48,
    'recall': 0.9230769230769231,
    'f1': 0.631578947368421,
    'accuracy': None}