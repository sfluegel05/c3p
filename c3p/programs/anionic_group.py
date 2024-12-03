"""
Classifies: CHEBI:64767 anionic group
"""
from rdkit import Chem

def is_anionic_group(smiles: str):
    """
    Determines if a molecule is an anionic group (a group that carries an overall negative charge).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anionic group, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of negative charges in the molecule
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() < 0:
            return True, f"Atom {atom.GetSymbol()} carries a negative charge"

    return False, "No atoms with negative charges found"

# Example usage
smiles_examples = [
    "C1=CC(=NC(N1[C@@H]2O[C@H](COP(O[C@H]3[C@H]([C@H](*)O[C@@H]3COP(*)(=O)[O-])O)([O-])=O)[C@H]([C@H]2O)OP(OC[C@H]4O[C@@H](N5C=CC(=NC5=O)N)[C@@H]([C@@H]4OP(OC[C@H]6O[C@@H](N7C8=C(N=C7)C(=NC=N8)N)[C@@H]([C@@H]6OP(OC[C@H]9O[C@@H](N%10C=CC(=NC%10=O)N)[C@@H]([C@@H]9OP(OC[C@H]%11O[C@@H](N%12C=CC(=NC%12=O)N)[C@@H]([C@@H]%11OP(OC[C@H]%13O[C@@H](N%14C%15=C(N=C%14)C(=NC=N%15)N)[C@@H]([C@@H]%13O)O)([O-])=O)O)([O-])=O)O)([O-])=O)O)([O-])=O)O)([O-])=O)=O)N",
    "C(C[C@@H](C[C@@H](CC(/C=C/C=C/C=C/CCCCC)O)O)O)(SCCNC(CCNC(=O)[C@@H](C(COP(OC[C@@H](C(*)=O)N*)(=O)[O-])(C)C)O)=O)=O",
    "N(C(C(NCC([O-])=O)=O)*)*",
    "C1(=O)NC(=NC2=C1N=CN2[C@@H]3O[C@H](COP(=O)(*)[O-])[C@@H](O*)[C@H]3O)NC(C([H])=O)O",
    "C(*)(=O)[C@@H](N*)CCC(=O)[O-]"
]

for smiles in smiles_examples:
    result, reason = is_anionic_group(smiles)
    print(f"SMILES: {smiles}\nIs anionic group: {result}\nReason: {reason}\n")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64767',
                          'name': 'anionic group',
                          'definition': 'A group that carries an overall '
                                        'negative charge.',
                          'parents': ['CHEBI:24433']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 9-10: malformed \\N character escape (<string>, line 1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}