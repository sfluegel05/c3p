"""
Classifies: CHEBI:25961 phenanthrenes
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_phenanthrenes(smiles: str):
    """
    Determines if a molecule is a phenanthrene or its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenanthrene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least three 6-membered rings
    if len([ring for ring in rings.AtomRings() if len(ring) == 6]) < 3:
        return False, "Less than three 6-membered rings found"

    # Find all aromatic 6-membered rings
    aromatic_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() for atom in atoms):
                aromatic_rings.append(ring)

    if len(aromatic_rings) < 3:
        return False, "Less than three aromatic 6-membered rings found"

    # Check for the phenanthrene structure (three fused benzene rings in a specific arrangement)
    phenanthrene_pattern = Chem.MolFromSmarts('c1ccc2c(c1)ccc3ccccc23')
    if not mol.HasSubstructMatch(phenanthrene_pattern):
        return False, "Phenanthrene core structure not found"

    return True, "Phenanthrene core structure found"

# Example usage
smiles = "O1[C@H]2C=Cc3c(ccc4ccccc34)[C@@H]12"
result, reason = is_phenanthrenes(smiles)
print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25961',
                          'name': 'phenanthrenes',
                          'definition': 'Any benzenoid aromatic compound that '
                                        'consists of a phenanthrene skeleton '
                                        'and its substituted derivatives '
                                        'thereof.',
                          'parents': ['CHEBI:33836', 'CHEBI:38032']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'False Less than three aromatic 6-membered rings found\n',
    'num_true_positives': 11,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 20,
    'precision': 1.0,
    'recall': 0.3548387096774194,
    'f1': 0.5238095238095238,
    'accuracy': None}