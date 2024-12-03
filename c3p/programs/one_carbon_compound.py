"""
Classifies: CHEBI:64708 one-carbon compound
"""
from rdkit import Chem

def is_one_carbon_compound(smiles: str):
    """
    Determines if a molecule is a one-carbon compound (contains a single carbon atom).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a one-carbon compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms in the molecule
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')

    if carbon_count == 1:
        return True, "Molecule contains a single carbon atom"
    else:
        return False, f"Molecule contains {carbon_count} carbon atoms"

# Examples for testing
print(is_one_carbon_compound("FC(F)(F)F"))  # tetrafluoromethane
print(is_one_carbon_compound("[H]C([H])([H])Br"))  # bromomethane
print(is_one_carbon_compound("[Na+].[Na+].OP([O-])(=O)C(Cl)(Cl)P(O)([O-])=O"))  # clodronic acid disodium salt
print(is_one_carbon_compound("[C-]#[O+]"))  # carbon monoxide
print(is_one_carbon_compound("[H]C([H])([H])[H]"))  # methane
print(is_one_carbon_compound("CS(O)(=O)=O"))  # methanesulfonic acid
print(is_one_carbon_compound("[H]C([H])([H])Cl"))  # chloromethane
print(is_one_carbon_compound("NCP(O)(O)=O"))  # (aminomethyl)phosphonic acid
print(is_one_carbon_compound("NC(=O)OP(O)(O)=O"))  # carbamoyl phosphate
print(is_one_carbon_compound("NC(O)=N"))  # carbamimidic acid
print(is_one_carbon_compound("[H]C([H])=O"))  # formaldehyde
print(is_one_carbon_compound("S=N#C"))  # thiofulminic acid


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64708',
                          'name': 'one-carbon compound',
                          'definition': 'An organic molecular entity '
                                        'containing a single carbon atom (C1).',
                          'parents': ['CHEBI:50860']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[02:21:03] Explicit valence for atom # 1 N, 5, is greater than '
             'permitted\n'
             '[02:21:03] Explicit valence for atom # 1 N, 5, is greater than '
             'permitted\n',
    'stdout': "(True, 'Molecule contains a single carbon atom')\n"
              "(True, 'Molecule contains a single carbon atom')\n"
              "(True, 'Molecule contains a single carbon atom')\n"
              "(True, 'Molecule contains a single carbon atom')\n"
              "(True, 'Molecule contains a single carbon atom')\n"
              "(True, 'Molecule contains a single carbon atom')\n"
              "(True, 'Molecule contains a single carbon atom')\n"
              "(True, 'Molecule contains a single carbon atom')\n"
              "(True, 'Molecule contains a single carbon atom')\n"
              "(True, 'Molecule contains a single carbon atom')\n"
              "(True, 'Molecule contains a single carbon atom')\n"
              "(False, 'Invalid SMILES string')\n",
    'num_true_positives': 11,
    'num_false_positives': 0,
    'num_true_negatives': 12,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9166666666666666,
    'f1': 0.9565217391304348,
    'accuracy': None}