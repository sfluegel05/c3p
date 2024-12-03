"""
Classifies: CHEBI:51307 diester
"""
from rdkit import Chem

def is_diester(smiles: str):
    """
    Determines if a molecule is a diester (contains two ester groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the ester functional group pattern
    ester_pattern = Chem.MolFromSmarts('C(=O)O')

    # Find all ester groups in the molecule
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if len(ester_matches) >= 2:
        return True, f"Contains {len(ester_matches)} ester groups"
    else:
        return False, f"Contains {len(ester_matches)} ester groups"

# Test cases
print(is_diester("CCC(=O)OC1=C(N=CC=C1OC)C(=O)N[C@@H](C)C(=O)O[C@@H](C)[C@@H](C)C1=C(C)C=CC=C1"))  # metarylpicoxamid
print(is_diester("COC(=O)c1c(Cl)c(Cl)c(C(=O)OC)c(Cl)c1Cl"))  # Dacthal
print(is_diester("[Cl-].[Cl-].COc1cc(C[C@H]2c3c(CC[N@+]2(C)CCCOC(=O)CCC(=O)OCCC[N@+]2(C)CCc4cc(OC)c(OC)c(OC)c4[C@H]2Cc2cc(OC)c(OC)c(OC)c2)cc(OC)c(OC)c3OC)cc(OC)c1OC"))  # meso-doxacurium chloride
print(is_diester("CC(C)OC(=O)c1cc(cc(c1)[N+]([O-])=O)C(=O)OC(C)C"))  # Nitrothal-isopropyl
print(is_diester("C\\C(\\C=C\\C=C(/C)C(=O)O[C@@H]1O[C@H](CO[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O)=C/C=C/C=C(\\C)/C=C/C=C(\\C)C(=O)O[C@@H]1O[C@H](CO[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O"))  # crocin-1
print(is_diester("[Cl-].[Cl-].COc1cc(C[C@H]2c3c(CC[N@+]2(C)CCCOC(=O)CCC(=O)OCCC[N@@+]2(C)CCc4cc(OC)c(OC)c(OC)c4[C@@H]2Cc2cc(OC)c(OC)c(OC)c2)cc(OC)c(OC)c3OC)cc(OC)c1OC"))  # (1S,2R,1'S,2'R)-doxacurium chloride
print(is_diester("CC[C@H](C)[C@H](NC(=O)[C@H]1CCCCN1C)C(=O)N(COC(=O)CC(C)C)[C@H](C[C@@H](OC(C)=O)c1nc(cs1)C(=O)N[C@H](C[C@H](C)C(O)=O)Cc1ccc(O)cc1)C(C)C"))  # Tubulysin A
print(is_diester("CC(C)CCCCCCOC(=O)c1ccccc1C(=O)OCCCCCCC(C)C"))  # diisononyl phthalate
print(is_diester("COC(=O)CCCCCCCC(=O)OC"))  # dimethyl azelate
print(is_diester("CCC1=C(c2ccc(OC(C)=O)cc2)C(C)(C)Oc2cc(OC(C)=O)ccc12"))  # 2,2-dimethyl-3-[4-(acetyloxy)phenyl]-4-ethyl-2H-1-benzopyran-7-ol acetate
print(is_diester("CC(=O)OCCCOC(C)=O"))  # 1,3-diacetoxypropane
print(is_diester("CC[C@H](C)[C@H](NC(=O)[C@H]1CCCCN1C)C(=O)N(COC(=O)CC(C)C)[C@H](C[C@@H](OC(C)=O)c1nc(cs1)C(=O)N[C@H](C[C@H](C)C(O)=O)Cc1ccccc1)C(C)C"))  # Tubulysin D


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51307',
                          'name': 'diester',
                          'definition': 'A diester is a compound containing '
                                        'two ester groups.',
                          'parents': ['CHEBI:35701']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Contains 2 ester groups')\n"
              "(True, 'Contains 2 ester groups')\n"
              "(True, 'Contains 2 ester groups')\n"
              "(True, 'Contains 2 ester groups')\n"
              "(True, 'Contains 2 ester groups')\n"
              "(True, 'Contains 2 ester groups')\n"
              "(True, 'Contains 3 ester groups')\n"
              "(True, 'Contains 2 ester groups')\n"
              "(True, 'Contains 2 ester groups')\n"
              "(True, 'Contains 2 ester groups')\n"
              "(True, 'Contains 2 ester groups')\n"
              "(True, 'Contains 3 ester groups')\n",
    'num_true_positives': 12,
    'num_false_positives': 1,
    'num_true_negatives': 11,
    'num_false_negatives': 0,
    'precision': 0.9230769230769231,
    'recall': 1.0,
    'f1': 0.9600000000000001,
    'accuracy': None}