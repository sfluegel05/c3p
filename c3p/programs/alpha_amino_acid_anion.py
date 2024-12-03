"""
Classifies: CHEBI:33558 alpha-amino-acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid_anion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino-acid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of carboxylate group (C(=O)[O-])
    carboxylate = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 0:
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if 'O' in neighbors and '[O-]' in neighbors:
                carboxylate = True
                break

    if not carboxylate:
        return False, "No carboxylate group found"

    # Check for the presence of alpha-amino group (N attached to the alpha carbon)
    alpha_amino = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    # Check if the carbon is alpha carbon (attached to carboxylate group)
                    carboxylate_neighbors = [n.GetSymbol() for n in neighbor.GetNeighbors()]
                    if 'C' in carboxylate_neighbors and 'O' in carboxylate_neighbors and '[O-]' in carboxylate_neighbors:
                        alpha_amino = True
                        break
            if alpha_amino:
                break

    if not alpha_amino:
        return False, "No alpha-amino group found"

    return True, "Alpha-amino-acid anion"

# Example usage
smiles_list = [
    "N[C@@H](CC(=O)OP([O-])([O-])=O)C([O-])=O",
    "NC1(CC1)C([O-])=O",
    "O=C([O-])[C@@H](N(O)O)CCCCCCCCSC",
    "C[C@@H]([NH2+]CCS([O-])(=O)=O)C([O-])=O",
    "N[C@@H](CC#N)C([O-])=O",
    "N[C@H](C([O-])=O)C(=O)COP([O-])([O-])=O",
    "CC(=O)N[C@@H](CS[*])C([O-])=O",
    "OC1CNC(C1)C([O-])=O",
    "N[C@@H](C[Se-])C([O-])=O",
    "C([C@@H](C([O-])=O)[NH3+])CCC[NH2+]CC([C@@H]([C@@H](CO)O)OP([O-])([O-])=O)=O",
    "CC(C)=CCNc1ncn(CC[C@H](N)C([O-])=O)c2ncnc12",
    "C[C@@H]([NH2+]CC([O-])=O)C([O-])=O",
    "NC(CC(N)=O)C([O-])=O",
    "[H]C(=O)c1ccc2oc3cc(=O)c(N)c(SC[C@H](NC(C)=O)C([O-])=O)c3nc2c1",
    "[O-]C(=O)[C@@H]1CCCN1",
    "NC(C[S-])C([O-])=O",
    "N[C@@H](C[SeH])C([O-])=O",
    "[NH3+][C@@H](CC([O-])=O)C([O-])=O",
    "O=C([O-])[C@@H](NO)CCCCCCCCSC",
    "[NH3+][C@@H](CCCC[NH2+][C@@H](CCC([O-])=O)C([O-])=O)C([O-])=O",
    "[C@H](C(=O)[O-])(N)*",
    "CSCCCC(N)C([O-])=O",
    "[O-]C(=O)C1CC(=O)CN1",
    "[NH3+][C@@H](C\\C(C=O)=C\\C=C(\\O)C([O-])=O)C([O-])=O",
    "CCCC[C@H](N)C([O-])=O"
]

for smiles in smiles_list:
    result, reason = is_alpha_amino_acid_anion(smiles)
    print(f"SMILES: {smiles} -> {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33558',
                          'name': 'alpha-amino-acid anion',
                          'definition': 'An amino-acid anion obtained by '
                                        'deprotonation of any alpha-amino '
                                        'acid.',
                          'parents': ['CHEBI:37022']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 8-9: malformed \\N character escape (<string>, line 1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}