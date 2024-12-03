"""
Classifies: CHEBI:23931 epoxy monocarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_epoxy_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is an epoxy monocarboxylic acid (Monocarboxylic acids containing at least one epoxy group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxy monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (COOH)
    carboxylic_acid = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('O') == 2 and neighbors.count('H') == 1:
                carboxylic_acid = True
                break

    if not carboxylic_acid:
        return False, "No carboxylic acid group found"

    # Check for epoxy group (three-membered ring with an oxygen atom)
    epoxy_group = False
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) == 3:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if any(atom.GetSymbol() == 'O' for atom in atoms):
                epoxy_group = True
                break

    if not epoxy_group:
        return False, "No epoxy group found"

    return True, "Epoxy monocarboxylic acid"

# Example usage
smiles = "CCCCCC1OC1C\C=C/CCCCCCCC(O)=O"  # Vernolic acid
result, reason = is_epoxy_monocarboxylic_acid(smiles)
print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23931',
                          'name': 'epoxy monocarboxylic acid',
                          'definition': 'Monocarboxylic acids containing at '
                                        'least one epoxy group.',
                          'parents': ['CHEBI:25384', 'CHEBI:32955']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 17-18: malformed \\N character escape (<string>, line '
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