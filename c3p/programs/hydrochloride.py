"""
Classifies: CHEBI:36807 hydrochloride
"""
from rdkit import Chem

def is_hydrochloride(smiles: str):
    """
    Determines if a molecule is a hydrochloride (a salt resulting from the reaction of hydrochloric acid with an organic base).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydrochloride, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of chloride ion
    chloride_present = any(atom.GetSymbol() == 'Cl' for atom in mol.GetAtoms())
    if not chloride_present:
        return False, "No chloride ion present"

    # Check for the presence of a positively charged organic base or a neutral amine
    positive_charge_present = any(atom.GetFormalCharge() == 1 for atom in mol.GetAtoms())
    neutral_amine_present = any(atom.GetSymbol() == 'N' and atom.GetFormalCharge() == 0 for atom in mol.GetAtoms())
    
    if not positive_charge_present and not neutral_amine_present:
        return False, "No positively charged organic base or neutral amine present"

    return True, "Molecule is a hydrochloride"

# Example usage:
smiles_examples = [
    "[H+].[Cl-].[H][C@@]12Cc3c(ccc(O)c3C(=O)C1=C(O)[C@]1(O)C(=O)C(C(N)=O)=C(O)[C@@H](N(C)C)[C@]1([H])C2)N(C)C",
    "N1=CC2=C(C(=C1)C)C(=CC=C2)S(N3[C@H](CNCCC3)C)(=O)=O.Cl.Cl",
    "O(C=1C(CCNCC=2C(O)=CC=CC2)=CC(OC)=C(C1)C#N)C.Cl",
    "Cl.Nc1c2ccccc2nc2ccccc12",
    "[H+].[H+].[Cl-].[Cl-].Cc1cc(nnc1NCCN1CCOCC1)-c1ccccc1",
    "Cl.CN(Cc1ccc(cc1)C(C)(C)C)Cc1cccc2ccccc12",
    "C=1(OC[C@H]2[C@](C=3C=CC(=CC3)F)(CCNC2)[H])C=C4C(OCO4)=CC1.[H]Cl",
    "O.OC[C@H](NC1=NC=2N(C=NC2C(=N1)NC3=CC(C4=NC=CC=C4)=CC=C3)C(C)C)CC.Cl",
    "Cl.[H][C@@]12[C@H](OC(=O)CCN(C)C)[C@H](OC(C)=O)[C@@]3(C)O[C@](C)(CC(=O)[C@]3(O)[C@@]1(C)[C@@H](O)CCC2(C)C)C=C",
    "Cl[H].Cl[H].NCCN"
]

for smiles in smiles_examples:
    result, reason = is_hydrochloride(smiles)
    print(f"SMILES: {smiles} -> {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36807',
                          'name': 'hydrochloride',
                          'definition': 'A salt formally resulting from the '
                                        'reaction of hydrochloric acid with an '
                                        'organic base.',
                          'parents': ['CHEBI:36094']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '[23:58:32] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[23:58:32] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[23:58:32] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[23:58:32] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[23:58:32] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[23:58:32] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[23:58:32] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[23:58:32] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[23:58:32] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[23:58:32] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[23:58:32] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[23:58:32] SMILES Parse Error: syntax error while parsing: '
             'C/C=C\x01/CN2CC[C@@]34C5=CC=CC=C5N6[C@@]([C@]1(C[C@]2([C@@]36OC7=C4C=C8C(=C7O)C(=O)O[C@@]9%10[C@@]%11(C[C@]%12(/C(=C\\C)/CN%11CC[C@@]8%10C%13=CC=CC=C%13N9[C@@]%12(C(=O)OC)[H])[H])[H])[H])[H])(C(=O)OC)[H]\n'
             '[23:58:32] SMILES Parse Error: Failed parsing SMILES '
             "'C/C=C\x01/CN2CC[C@@]34C5=CC=CC=C5N6[C@@]([C@]1(C[C@]2([C@@]36OC7=C4C=C8C(=C7O)C(=O)O[C@@]9%10[C@@]%11(C[C@]%12(/C(=C\\C)/CN%11CC[C@@]8%10C%13=CC=CC=C%13N9[C@@]%12(C(=O)OC)[H])[H])[H])[H])[H])(C(=O)OC)[H]' "
             'for input: '
             "'C/C=C\x01/CN2CC[C@@]34C5=CC=CC=C5N6[C@@]([C@]1(C[C@]2([C@@]36OC7=C4C=C8C(=C7O)C(=O)O[C@@]9%10[C@@]%11(C[C@]%12(/C(=C\\C)/CN%11CC[C@@]8%10C%13=CC=CC=C%13N9[C@@]%12(C(=O)OC)[H])[H])[H])[H])[H])(C(=O)OC)[H]'\n"
             '[23:58:32] Explicit valence for atom # 1 Al, 9, is greater than '
             'permitted\n',
    'stdout': 'SMILES: '
              '[H+].[Cl-].[H][C@@]12Cc3c(ccc(O)c3C(=O)C1=C(O)[C@]1(O)C(=O)C(C(N)=O)=C(O)[C@@H](N(C)C)[C@]1([H])C2)N(C)C '
              '-> True, Reason: Molecule is a hydrochloride\n'
              'SMILES: '
              'N1=CC2=C(C(=C1)C)C(=CC=C2)S(N3[C@H](CNCCC3)C)(=O)=O.Cl.Cl -> '
              'True, Reason: Molecule is a hydrochloride\n'
              'SMILES: O(C=1C(CCNCC=2C(O)=CC=CC2)=CC(OC)=C(C1)C#N)C.Cl -> '
              'True, Reason: Molecule is a hydrochloride\n'
              'SMILES: Cl.Nc1c2ccccc2nc2ccccc12 -> True, Reason: Molecule is a '
              'hydrochloride\n'
              'SMILES: [H+].[H+].[Cl-].[Cl-].Cc1cc(nnc1NCCN1CCOCC1)-c1ccccc1 '
              '-> True, Reason: Molecule is a hydrochloride\n'
              'SMILES: Cl.CN(Cc1ccc(cc1)C(C)(C)C)Cc1cccc2ccccc12 -> True, '
              'Reason: Molecule is a hydrochloride\n'
              'SMILES: '
              'C=1(OC[C@H]2[C@](C=3C=CC(=CC3)F)(CCNC2)[H])C=C4C(OCO4)=CC1.[H]Cl '
              '-> True, Reason: Molecule is a hydrochloride\n'
              'SMILES: '
              'O.OC[C@H](NC1=NC=2N(C=NC2C(=N1)NC3=CC(C4=NC=CC=C4)=CC=C3)C(C)C)CC.Cl '
              '-> True, Reason: Molecule is a hydrochloride\n'
              'SMILES: '
              'Cl.[H][C@@]12[C@H](OC(=O)CCN(C)C)[C@H](OC(C)=O)[C@@]3(C)O[C@](C)(CC(=O)[C@]3(O)[C@@]1(C)[C@@H](O)CCC2(C)C)C=C '
              '-> True, Reason: Molecule is a hydrochloride\n'
              'SMILES: Cl[H].Cl[H].NCCN -> True, Reason: Molecule is a '
              'hydrochloride\n',
    'num_true_positives': 38,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}