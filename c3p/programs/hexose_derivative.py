"""
Classifies: CHEBI:63385 hexose derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_hexose_derivative(smiles: str):
    """
    Determines if a molecule is a hexose derivative (a monosaccharide derivative that is formally obtained from a hexose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define hexose backbones for pyranose and furanose forms
    hexopyranose = Chem.MolFromSmarts("C1(C(C(C(C(O1)O)O)O)O)O")
    hexofuranose = Chem.MolFromSmarts("C1(C(C(O1)O)O)O")

    if mol.HasSubstructMatch(hexopyranose) or mol.HasSubstructMatch(hexofuranose):
        return True, "Molecule is a hexose derivative"
    
    # Check for derivatives by looking for common functional groups attached to a hexose backbone
    phosphate = Chem.MolFromSmarts("OP(O)(O)=O")
    if mol.HasSubstructMatch(phosphate):
        if mol.HasSubstructMatch(hexopyranose) or mol.HasSubstructMatch(hexofuranose):
            return True, "Molecule is a hexose derivative with phosphate group"
    
    # Check for other common hexose derivatives
    other_derivatives = [
        Chem.MolFromSmarts("C1(C(C(C(C(O1)O)O)O)O)O"),  # Hexopyranose
        Chem.MolFromSmarts("C1(C(C(O1)O)O)O"),  # Hexofuranose
        Chem.MolFromSmarts("C(C(C(C(C(C=O)O)O)O)O)O"),  # Aldohexose
        Chem.MolFromSmarts("C(C(C(C(C(C=O)O)O)O)O)O"),  # Ketohexose
    ]
    
    for derivative in other_derivatives:
        if mol.HasSubstructMatch(derivative):
            return True, "Molecule is a hexose derivative"
    
    return False, "No hexose backbone found"

# Test cases
print(is_hexose_derivative("[H][C@@]1(O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)Oc1c[nH]c2ccc(Br)c(Cl)c12"))  # 5-bromo-4-chloro-3-indolyl alpha-D-glucoside
print(is_hexose_derivative("OC[C@@H]1O[C@@H](OP(O)(O)=O)[C@@H](O)[C@H](O)[C@@H]1O"))  # alpha-L-galactose 1-phosphate
print(is_hexose_derivative("O[C@H](COP(O)(O)=O)[C@H](O)[C@@H](O)C(=O)COP(O)(O)=O"))  # D-sorbose 1,6-bisphosphate
print(is_hexose_derivative("CCCCCCCC(C)CCCC(C)CCCC(C)CCCC(C)CCCC(C)CCCOP(O)(=O)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O"))  # beta-D-mannosyl 4,8,12,16,20-pentamethylheptacosyl phosphate
print(is_hexose_derivative("OCC1(OP(O)(O)=O)O[C@H](COP(O)(O)=O)[C@@H](O)[C@@H]1O"))  # D-fructofuranose 2,6-bisphosphate
print(is_hexose_derivative("OC1C[C@@H](O)[C@H](O)[C@@H](COP(O)(O)=O)O1"))  # 2-deoxy-D-glucopyranose 6-phosphate
print(is_hexose_derivative("[H]C(=O)[C@H](OP(O)(O)=O)[C@@H](O)[C@H](O)[C@H](O)COP(O)(O)=O"))  # D-glucose 1,6-bisphosphate
print(is_hexose_derivative("O[C@H](COP(O)(O)=O)[C@@H](O)[C@H](O)C(=O)COP(O)(O)=O"))  # keto-D-fructose 1,6-bisphosphate
print(is_hexose_derivative("OC1O[C@H](COP(O)(O)=O)C(O)C(O)C1O"))  # D-hexopyranose 6-phosphate
print(is_hexose_derivative("O[C@H]1[C@H](O)C(O)(COP(O)(O)=O)O[C@@H]1COP(O)(O)=O"))  # D-fructofuranose 1,6-bisphosphate


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:63385',
                          'name': 'hexose derivative',
                          'definition': 'A monosaccharide derivative that is '
                                        'formally obtained from a hexose.',
                          'parents': ['CHEBI:63367']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(False, 'No hexose backbone found')\n"
              "(False, 'No hexose backbone found')\n"
              "(False, 'No hexose backbone found')\n"
              "(False, 'No hexose backbone found')\n"
              "(False, 'No hexose backbone found')\n"
              "(False, 'No hexose backbone found')\n"
              "(True, 'Molecule is a hexose derivative')\n"
              "(False, 'No hexose backbone found')\n"
              "(False, 'No hexose backbone found')\n"
              "(False, 'No hexose backbone found')\n",
    'num_true_positives': 1,
    'num_false_positives': 1,
    'num_true_negatives': 9,
    'num_false_negatives': 9,
    'precision': 0.5,
    'recall': 0.1,
    'f1': 0.16666666666666669,
    'accuracy': None}