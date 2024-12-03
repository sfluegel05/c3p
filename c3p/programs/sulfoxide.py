"""
Classifies: CHEBI:22063 sulfoxide
"""
from rdkit import Chem

def is_sulfoxide(smiles: str):
    """
    Determines if a molecule is a sulfoxide (R2S=O or R2C=S=O where R is not H).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfoxide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of sulfoxide group (S=O)
    sulfoxide_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'S':
            oxygen_neighbors = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() == 2.0]
            if len(oxygen_neighbors) == 1:
                carbon_neighbors = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == 'C']
                if len(carbon_neighbors) >= 1:
                    sulfoxide_found = True
                    break

    if not sulfoxide_found:
        return False, "No sulfoxide group (S=O) found"

    # Check that sulfur has two organic substituents (R groups) and no hydrogen atoms
    r_count = 0
    for neighbor in atom.GetNeighbors():
        if neighbor.GetSymbol() != 'O' and neighbor.GetSymbol() != 'H':
            r_count += 1
        if neighbor.GetSymbol() == 'H':
            return False, "Sulfur atom has hydrogen atom as a neighbor"

    if r_count < 2:
        return False, "Sulfur atom does not have two organic substituents"

    return True, "Valid sulfoxide group found"

# Example usage
smiles_list = [
    "S(=O)(C(SSCCC)CC)/C=C/C",
    "S(=O)(C(SS/C=C/C)CC)/C=C\\C",
    "OC1CC(O)C(O)C(O)CS(=O)CC(O)C(O)C(O)C(O)C1",
    "CNC(=O)Oc1cc(C)c(c(C)c1)S(C)=O",
    "O1[C@@H]([C@H]([C@@H]([C@H]([C@@H]1SC(CCCCCCCS(C)=O)=NOS(O)(=O)=O)O)O)O)CO",
    "C=1C(=CC=CC1)[S@@+](C)[O-]",
    "COc1ccnc(CS(=O)c2nc3cc(OC(F)F)ccc3[nH]2)c1OC",
    "S(=O)(C(SSC)CC)CC=C",
    "C=1C=C(C=C2C1[C@]3(CC[C@@]4([C@H](CC[C@]4([C@@]3([C@@H](C2)CCCCCCCCCS(=O)CCCC(C(F)(F)F)(F)F)[H])[H])O)C)[H])O",
    "CC(C)OC(=O)C(C(=O)OC(C)C)=C1SCCS1=O",
    "N(=C=S)CCCCS(=O)C",
    "[S@](CC1=C(C)C(OCCCOC)=CC=N1)(=O)C2=NC3=C(N2)C=CC=C3",
    "CCCS(SCCC)=O",
    "COc1ccc2[nH]c(nc2c1)S(=O)Cc1ncc(C)c(OC)c1C",
    "C[N+](C)(C)[C@@H](Cc1c[nH]c([Se]C[C@H](NC(=O)CC[C@H](N)C(O)=O)C(O)=O)n1)C([O-])=O"
]

for smiles in smiles_list:
    result, reason = is_sulfoxide(smiles)
    print(f"SMILES: {smiles} -> {result}, {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22063',
                          'name': 'sulfoxide',
                          'definition': 'An organosulfur compound having the '
                                        'structure R2S=O or R2C=S=O (R =/= H).',
                          'parents': ['CHEBI:33261']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: S(=O)(C(SSCCC)CC)/C=C/C -> True, Valid sulfoxide group '
              'found\n'
              'SMILES: S(=O)(C(SS/C=C/C)CC)/C=C\\C -> True, Valid sulfoxide '
              'group found\n'
              'SMILES: OC1CC(O)C(O)C(O)CS(=O)CC(O)C(O)C(O)C(O)C1 -> True, '
              'Valid sulfoxide group found\n'
              'SMILES: CNC(=O)Oc1cc(C)c(c(C)c1)S(C)=O -> True, Valid sulfoxide '
              'group found\n'
              'SMILES: '
              'O1[C@@H]([C@H]([C@@H]([C@H]([C@@H]1SC(CCCCCCCS(C)=O)=NOS(O)(=O)=O)O)O)O)CO '
              '-> True, Valid sulfoxide group found\n'
              'SMILES: C=1C(=CC=CC1)[S@@+](C)[O-] -> False, No sulfoxide group '
              '(S=O) found\n'
              'SMILES: COc1ccnc(CS(=O)c2nc3cc(OC(F)F)ccc3[nH]2)c1OC -> True, '
              'Valid sulfoxide group found\n'
              'SMILES: S(=O)(C(SSC)CC)CC=C -> True, Valid sulfoxide group '
              'found\n'
              'SMILES: '
              'C=1C=C(C=C2C1[C@]3(CC[C@@]4([C@H](CC[C@]4([C@@]3([C@@H](C2)CCCCCCCCCS(=O)CCCC(C(F)(F)F)(F)F)[H])[H])O)C)[H])O '
              '-> True, Valid sulfoxide group found\n'
              'SMILES: CC(C)OC(=O)C(C(=O)OC(C)C)=C1SCCS1=O -> True, Valid '
              'sulfoxide group found\n'
              'SMILES: N(=C=S)CCCCS(=O)C -> True, Valid sulfoxide group found\n'
              'SMILES: [S@](CC1=C(C)C(OCCCOC)=CC=N1)(=O)C2=NC3=C(N2)C=CC=C3 -> '
              'True, Valid sulfoxide group found\n'
              'SMILES: CCCS(SCCC)=O -> True, Valid sulfoxide group found\n'
              'SMILES: COc1ccc2[nH]c(nc2c1)S(=O)Cc1ncc(C)c(OC)c1C -> True, '
              'Valid sulfoxide group found\n'
              'SMILES: '
              'C[N+](C)(C)[C@@H](Cc1c[nH]c([Se]C[C@H](NC(=O)CC[C@H](N)C(O)=O)C(O)=O)n1)C([O-])=O '
              '-> False, No sulfoxide group (S=O) found\n',
    'num_true_positives': 13,
    'num_false_positives': 0,
    'num_true_negatives': 15,
    'num_false_negatives': 2,
    'precision': 1.0,
    'recall': 0.8666666666666667,
    'f1': 0.9285714285714286,
    'accuracy': None}