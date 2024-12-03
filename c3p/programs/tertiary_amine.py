"""
Classifies: CHEBI:32876 tertiary amine
"""
from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for nitrogen atoms
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'N']
    if not nitrogen_atoms:
        return False, "No nitrogen atoms found"

    for atom in nitrogen_atoms:
        # Check if nitrogen is connected to exactly three carbon atoms
        neighbors = atom.GetNeighbors()
        if len(neighbors) == 3 and all(neighbor.GetSymbol() == 'C' for neighbor in neighbors):
            return True, "Tertiary amine found"
    
    return False, "No tertiary amine found"

# Example usage:
smiles_list = [
    "CC(CN1CCCCC1)Cc1ccc(cc1)C(C)(C)C",  # 1-[3-(4-tert-butylphenyl)-2-methylpropyl]piperidine
    "CCC#CC(CN(C)C)(C)C",  # N,N,2,2-tetramethyl-3-hexyn-1-amine
    "CC(C)N(C(C)C)C(=O)SCC(Cl)=C(Cl)Cl",  # Tri-allate
    "CCN(CC)CCOC(=O)C1(CCCCC1)C1CCCCC1",  # dicyclomine
    "[H]C(=C([H])c1cc[n+](C)cc1)c1ccc(cc1)N(CCCCCCCCCCCCCCCC)CCCCCCCCCCCCCCCC",  # 4-(4-dihexadecylaminostyryl)-N-methylpyridium
    "CCSC(=O)N(CC)CC",  # Ethiolate
    "CN1[C@H]2C[C@@H](C[C@@H]1[C@@H](O)C2)OC(=O)C(CO)c1ccccc1",  # (6S)-6-hydroxyhyoscyamine
    "CCCSC(=O)N(CCC)CCC",  # S-propyl dipropylcarbamothioate
    "CCN(CCCc1ccccc1)CCCc1ccccc1",  # alverine
    "C1=CC=C(C=C1)N(C2=CC=CC=C2)C3=CC=C(C=C3)C=NNC(=O)CCCCCCC(=O)NO"  # LSM-6269
]

for smiles in smiles_list:
    result, reason = is_tertiary_amine(smiles)
    print(f"SMILES: {smiles} -> {result}, {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:32876',
                          'name': 'tertiary amine',
                          'definition': 'A compound formally derived from '
                                        'ammonia by replacing three hydrogen '
                                        'atoms by hydrocarbyl groups.',
                          'parents': ['CHEBI:32952', 'CHEBI:50996']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: CC(CN1CCCCC1)Cc1ccc(cc1)C(C)(C)C -> True, Tertiary '
              'amine found\n'
              'SMILES: CCC#CC(CN(C)C)(C)C -> True, Tertiary amine found\n'
              'SMILES: CC(C)N(C(C)C)C(=O)SCC(Cl)=C(Cl)Cl -> True, Tertiary '
              'amine found\n'
              'SMILES: CCN(CC)CCOC(=O)C1(CCCCC1)C1CCCCC1 -> True, Tertiary '
              'amine found\n'
              'SMILES: '
              '[H]C(=C([H])c1cc[n+](C)cc1)c1ccc(cc1)N(CCCCCCCCCCCCCCCC)CCCCCCCCCCCCCCCC '
              '-> True, Tertiary amine found\n'
              'SMILES: CCSC(=O)N(CC)CC -> True, Tertiary amine found\n'
              'SMILES: '
              'CN1[C@H]2C[C@@H](C[C@@H]1[C@@H](O)C2)OC(=O)C(CO)c1ccccc1 -> '
              'True, Tertiary amine found\n'
              'SMILES: CCCSC(=O)N(CCC)CCC -> True, Tertiary amine found\n'
              'SMILES: CCN(CCCc1ccccc1)CCCc1ccccc1 -> True, Tertiary amine '
              'found\n'
              'SMILES: '
              'C1=CC=C(C=C1)N(C2=CC=CC=C2)C3=CC=C(C=C3)C=NNC(=O)CCCCCCC(=O)NO '
              '-> True, Tertiary amine found\n',
    'num_true_positives': 10,
    'num_false_positives': 10,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.5,
    'recall': 1.0,
    'f1': 0.6666666666666666,
    'accuracy': None}