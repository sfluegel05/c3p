"""
Classifies: CHEBI:231691 sialylglycoconjugate anion
"""
from rdkit import Chem

def is_sialylglycoconjugate_anion(smiles: str):
    """
    Determines if a molecule is a sialylglycoconjugate anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sialylglycoconjugate anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of sialic acid (N-acetylneuraminic acid) residues
    sialic_acid_smarts = Chem.MolFromSmarts('[C@@H]1([C@H](NC(C)=O)[C@H](O)[C@H](O1)CO)C(O)=O')
    if not mol.HasSubstructMatch(sialic_acid_smarts):
        return False, "No sialic acid residues found"

    # Check for the presence of an anionic carboxylate group
    carboxylate_smarts = Chem.MolFromSmarts('C(=O)[O-]')
    if not mol.HasSubstructMatch(carboxylate_smarts):
        return False, "No anionic carboxylate group found"

    # Check for the presence of glycosidic linkages
    glycosidic_linkage_smarts = Chem.MolFromSmarts('O[C@H]1[C@H](O)[C@@H](O)[C@@H](CO)O1')
    if not mol.HasSubstructMatch(glycosidic_linkage_smarts):
        return False, "No glycosidic linkages found"

    return True, "Molecule is a sialylglycoconjugate anion"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:231691',
                          'name': 'sialylglycoconjugate anion',
                          'definition': 'A sialylglycoconjugate where the '
                                        'composition of the glycoconjugate '
                                        'represented with an R is not defined, '
                                        'can be an oligosaccharide, a '
                                        'glycoprotein or a glycolipid.',
                          'parents': ['CHEBI:78616']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 29,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}