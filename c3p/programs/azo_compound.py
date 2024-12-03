"""
Classifies: CHEBI:37533 azo compound
"""
from rdkit import Chem

def is_azo_compound(smiles: str):
    """
    Determines if a molecule is an azo compound (R-N=N-R).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an azo compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the azo group (N=N)
    pattern = Chem.MolFromSmarts("N=N")
    if mol.HasSubstructMatch(pattern):
        return True, "Contains azo group (N=N)"
    else:
        return False, "Does not contain azo group (N=N)"

# Example usage:
smiles_list = [
    "C1(=CC(=C(C=C1OC)/N=N/C2=CC(=C(C=C2)OC)S(O)(=O)=O)C)N",
    "CN(C)CCC1=CC=CC=C1\\N=N\\C1=CC=C(O)C(O)=C1",
    "C1=C(C2=C(C(=C1)N=NC3=C4C5=C(C=C3)NC(NC5=CC=C4)(C)C)C=CC=C2)N=NC6=CC=CC=C6",
    "Oc1ccc(cc1\\N=N\\c1c(O)c2ccc(Nc3ncnc(Nc4ccc5c(O)c(\\N=N\\c6cc(ccc6O)S(O)(=O)=O)c(cc5c4S(O)(=O)=O)S(O)(=O)=O)n3)c(c2cc1S(O)(=O)=O)S(O)(=O)=O)S(O)(=O)=O",
    "N(S(=O)(=O)C1=CC=C(C=C1)/N=N/C2=CC(C(O)=O)=C(C=C2)O)C=3N=C(C)C=C(N3)C",
    "Nc1ccc(c(N)c1)\\N=N\\c1ccc(cc1)S(N)(=O)=O",
    "FC1=CC=CC(F)=C1N=NC1=NN(CC2=C(F)C=CC=C2C(F)(F)F)C=C1",
    "N(C=1C=2C(C(=CC1)/N=N/C=3C=4C(C=C(C3)S(=O)(=O)O)=CC(=CC4O)S(=O)(=O)O)=CC=CC2S(=O)(=O)O)C5=CC=CC=C5",
    "FC1=CC=CC(F)=C1\\N=N\\C1=NN(CC2=C(F)C=CC=C2C(F)(F)F)C=C1",
    "COP(=S)(OC)Oc1ccc(cc1)\\N=N\\c1ccc(Cl)cc1",
    "OC(=O)c1cc(ccc1O)N=Nc1ccc(O)c(c1)C(O)=O",
    "C1=CC2C(=CN=NC2=NN)C=C1",
    "Oc1ccc2ccccc2c1\\N=N\\c1ccccc1"
]

for smiles in smiles_list:
    result, reason = is_azo_compound(smiles)
    print(f"SMILES: {smiles} -> {result}, {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37533',
                          'name': 'azo compound',
                          'definition': 'Derivatives of diazene with the '
                                        "general structure R-N=N-R'.",
                          'parents': ['CHEBI:51143']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 18-19: malformed \\N character escape (<string>, line '
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