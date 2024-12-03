"""
Classifies: CHEBI:36326 octadecanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_octadecanoid(smiles: str):
    """
    Determines if a molecule is an octadecanoid (unsaturated C18 fatty acids and skeletally related compounds).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octadecanoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule has 18 carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if num_carbons != 18:
        return False, "Molecule does not have 18 carbon atoms"

    # Check for the presence of at least one double bond (unsaturation)
    has_double_bond = any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE for bond in mol.GetBonds())
    if not has_double_bond:
        return False, "Molecule does not have any double bonds (unsaturation)"

    # Check for the presence of carboxylic acid group (COOH)
    carboxylic_acid = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "Molecule does not have a carboxylic acid group (COOH)"

    return True, "Molecule is an octadecanoid"

# Example usage:
# smiles = "C1(C(C/C=C\CCCCCO)O1)CCCCCCCC(=O)O"
# print(is_octadecanoid(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36326',
                          'name': 'octadecanoid',
                          'definition': 'Unsaturated C18 fatty acids and '
                                        'skeletally related compounds.',
                          'parents': ['CHEBI:15904', 'CHEBI:27208']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 11,
    'num_false_positives': 2,
    'num_true_negatives': 12,
    'num_false_negatives': 3,
    'precision': 0.8461538461538461,
    'recall': 0.7857142857142857,
    'f1': 0.8148148148148148,
    'accuracy': None}