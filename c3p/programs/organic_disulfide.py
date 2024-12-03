"""
Classifies: CHEBI:35489 organic disulfide
"""
from rdkit import Chem

def is_organic_disulfide(smiles: str):
    """
    Determines if a molecule is an organic disulfide (compounds of structure RSSR in which R and R' are organic groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic disulfide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all sulfur atoms in the molecule
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'S']

    # Check for disulfide bonds (S-S)
    disulfide_bonds = []
    for i, atom1 in enumerate(sulfur_atoms):
        for atom2 in sulfur_atoms[i+1:]:
            if mol.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx()) is not None:
                disulfide_bonds.append((atom1, atom2))

    if not disulfide_bonds:
        return False, "No disulfide bonds (S-S) found"

    # Check if each sulfur in the disulfide bond is connected to an organic group
    for s1, s2 in disulfide_bonds:
        s1_neighbors = [n for n in s1.GetNeighbors() if n.GetSymbol() != 'H' and n.GetIdx() != s2.GetIdx()]
        s2_neighbors = [n for n in s2.GetNeighbors() if n.GetSymbol() != 'H' and n.GetIdx() != s1.GetIdx()]

        if not s1_neighbors or not s2_neighbors:
            return False, "Disulfide bond not connected to organic groups"

        if not all(neigh.GetSymbol() in ['C', 'N', 'O', 'P', 'S'] for neigh in s1_neighbors + s2_neighbors):
            return False, "Disulfide bond connected to non-organic groups"

    return True, "Organic disulfide detected"

# Example usage:
smiles_list = [
    "CSCSSC",  # Methyl (methylthio)methyl disulfide
    "S1S[C@]23N([C@@H]4[C@@H](OC(=O)CC5=CC=CC=C5)C=CC=C4C2)C([C@]16N([C@@H]7[C@@H](OC(=O)[C@H](O)C8=CC=CC=C8)C=COC=C7C6)C3=O)=O",  # Emethallicin A
    "O([C@@H]1[C@@H]([C@H]([C@H]([C@H](O1)C)N)O)NC(C)=O)CCCCCCSSCCCCCCO[C@@H]2[C@@H]([C@H]([C@H]([C@H](O2)C)N)O)NC(C)=O",  # 6,6'-dithiodi[1-(2-acetamido-4-amino-2,4-dideoxy-alpha-D-fucosyloxy)hexane]
    "CC(C)(CO)[C@@H](O)C(=O)NCCC(=O)NCCSSCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)CO",  # pantethine
    "C(CCCSSCCCC(O)=O)(O)=O",  # 4,4'-disulfanyldibutanoic acid
    "O[C@H]1CSSC[C@H]1O",  # cis-1,2-dithiane-4,5-diol
    "C(SSCc1ccccc1)c1ccccc1",  # dibenzyl disulfide
    "[H][C@]12C[C@@]34SS[C@]5(C[C@]6([H])[C@@H](O)CC[C@H](O)[C@@]6([H])N5C3=O)C(=O)N4[C@]1([H])[C@@H](O)CC[C@@H]2O",  # rostratin A
    "N#CSSC#N",  # thiocyanogen
    "S(SCC)C(SC)C",  # Ethyl 1-(methylthio)ethyl disulfide
    "OC(=O)C1=C2SSC2=CC=CC1=O",  # tropodithietic acid
    "[H][C@]12CN3C4=C([C@@H](COC(N)=O)[C@@]3(OC)[C@@]1([H])N2)C(=O)C(NCCSSC1=CC=C(C=C1)[N+]([O-])=O)=C(C)C4=O",  # BMY-25067
    "C(=O)([C@@H](N)CSSC[C@@H](C(=O)O)N)NC=1C=CC2=C(C1)C=CC=C2",  # L-cystine mono-2-naphthylamide
    "N[C@@H](CCC(=O)N[C@H]1CSSC[C@H](NC(=O)CC[C@H](N)C(O)=O)C(=O)NCC(=O)NCCCNCCCCNC(=O)CNC1=O)C(O)=O"  # trypanothione disulfide
]

for smiles in smiles_list:
    result, reason = is_organic_disulfide(smiles)
    print(f"SMILES: {smiles}\nResult: {result}\nReason: {reason}\n")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35489',
                          'name': 'organic disulfide',
                          'definition': 'Compounds of structure RSSR in which '
                                        "R and R' are organic groups.",
                          'parents': ['CHEBI:33261', 'CHEBI:48343']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: CSCSSC\n'
              'Result: True\n'
              'Reason: Organic disulfide detected\n'
              '\n'
              'SMILES: '
              'S1S[C@]23N([C@@H]4[C@@H](OC(=O)CC5=CC=CC=C5)C=CC=C4C2)C([C@]16N([C@@H]7[C@@H](OC(=O)[C@H](O)C8=CC=CC=C8)C=COC=C7C6)C3=O)=O\n'
              'Result: True\n'
              'Reason: Organic disulfide detected\n'
              '\n'
              'SMILES: '
              'O([C@@H]1[C@@H]([C@H]([C@H]([C@H](O1)C)N)O)NC(C)=O)CCCCCCSSCCCCCCO[C@@H]2[C@@H]([C@H]([C@H]([C@H](O2)C)N)O)NC(C)=O\n'
              'Result: True\n'
              'Reason: Organic disulfide detected\n'
              '\n'
              'SMILES: '
              'CC(C)(CO)[C@@H](O)C(=O)NCCC(=O)NCCSSCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)CO\n'
              'Result: True\n'
              'Reason: Organic disulfide detected\n'
              '\n'
              'SMILES: C(CCCSSCCCC(O)=O)(O)=O\n'
              'Result: True\n'
              'Reason: Organic disulfide detected\n'
              '\n'
              'SMILES: O[C@H]1CSSC[C@H]1O\n'
              'Result: True\n'
              'Reason: Organic disulfide detected\n'
              '\n'
              'SMILES: C(SSCc1ccccc1)c1ccccc1\n'
              'Result: True\n'
              'Reason: Organic disulfide detected\n'
              '\n'
              'SMILES: '
              '[H][C@]12C[C@@]34SS[C@]5(C[C@]6([H])[C@@H](O)CC[C@H](O)[C@@]6([H])N5C3=O)C(=O)N4[C@]1([H])[C@@H](O)CC[C@@H]2O\n'
              'Result: True\n'
              'Reason: Organic disulfide detected\n'
              '\n'
              'SMILES: N#CSSC#N\n'
              'Result: True\n'
              'Reason: Organic disulfide detected\n'
              '\n'
              'SMILES: S(SCC)C(SC)C\n'
              'Result: True\n'
              'Reason: Organic disulfide detected\n'
              '\n'
              'SMILES: OC(=O)C1=C2SSC2=CC=CC1=O\n'
              'Result: True\n'
              'Reason: Organic disulfide detected\n'
              '\n'
              'SMILES: '
              '[H][C@]12CN3C4=C([C@@H](COC(N)=O)[C@@]3(OC)[C@@]1([H])N2)C(=O)C(NCCSSC1=CC=C(C=C1)[N+]([O-])=O)=C(C)C4=O\n'
              'Result: True\n'
              'Reason: Organic disulfide detected\n'
              '\n'
              'SMILES: '
              'C(=O)([C@@H](N)CSSC[C@@H](C(=O)O)N)NC=1C=CC2=C(C1)C=CC=C2\n'
              'Result: True\n'
              'Reason: Organic disulfide detected\n'
              '\n'
              'SMILES: '
              'N[C@@H](CCC(=O)N[C@H]1CSSC[C@H](NC(=O)CC[C@H](N)C(O)=O)C(=O)NCC(=O)NCCCNCCCCNC(=O)CNC1=O)C(O)=O\n'
              'Result: True\n'
              'Reason: Organic disulfide detected\n'
              '\n',
    'num_true_positives': 14,
    'num_false_positives': 0,
    'num_true_negatives': 14,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}