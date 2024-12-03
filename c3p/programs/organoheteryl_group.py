"""
Classifies: CHEBI:33456 organoheteryl group
"""
from rdkit import Chem

def is_organoheteryl_group(smiles: str):
    """
    Determines if a molecule is an organoheteryl group (a univalent group containing carbon which has its free valence at an atom other than carbon).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoheteryl group, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C']
    if not carbon_atoms:
        return False, "No carbon atoms found"
    
    # Check for univalent group (free valence)
    free_valence_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() != 'C' and atom.GetNumImplicitHs() > 0]
    if not free_valence_atoms:
        return False, "No free valence at an atom other than carbon"
    
    # Check if free valence is at an atom other than carbon
    for atom in free_valence_atoms:
        if any(neighbor.GetSymbol() == 'C' for neighbor in atom.GetNeighbors()):
            return True, "Organoheteryl group detected"
    
    return False, "Free valence is not at an atom other than carbon"

# Examples for testing
print(is_organoheteryl_group("OC(=O)[C@@H](N*)CC=1N=CNC1"))  # N(2)-L-histidino group
print(is_organoheteryl_group("O=C(N)[C@@H](N*)CCSC"))  # L-methionine amide residue
print(is_organoheteryl_group("O=C(O)C(N*)CCC(=O)N"))  # N(2)-glutamino group
print(is_organoheteryl_group("CO*"))  # methoxy group
print(is_organoheteryl_group("C(N)(=O)[C@@H](N*)CC=1N=CNC1"))  # L-histidine amide residue
print(is_organoheteryl_group("O=C(O)C(N*)CC=1C=CC(=CC1)O"))  # tyrosino group
print(is_organoheteryl_group("N(C(C(N)=O)*)*"))  # C-terminal alpha-amino-acid amide residue
print(is_organoheteryl_group("C1(=CNC2=C1C=CC=C2)C[C@H](C(=O)O)N*"))  # D-tryptophano group
print(is_organoheteryl_group("NC(=O)N-*"))  # carbamoylamino group
print(is_organoheteryl_group("O=C(O*)C(=O)O"))  # oxalooxy group


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33456',
                          'name': 'organoheteryl group',
                          'definition': 'A univalent group containing carbon '
                                        'which has its free valence at an atom '
                                        'other than carbon.',
                          'parents': ['CHEBI:51447']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Organoheteryl group detected')\n"
              "(True, 'Organoheteryl group detected')\n"
              "(True, 'Organoheteryl group detected')\n"
              "(False, 'No free valence at an atom other than carbon')\n"
              "(True, 'Organoheteryl group detected')\n"
              "(True, 'Organoheteryl group detected')\n"
              "(True, 'Organoheteryl group detected')\n"
              "(True, 'Organoheteryl group detected')\n"
              "(True, 'Organoheteryl group detected')\n"
              "(True, 'Organoheteryl group detected')\n",
    'num_true_positives': 9,
    'num_false_positives': 10,
    'num_true_negatives': 0,
    'num_false_negatives': 1,
    'precision': 0.47368421052631576,
    'recall': 0.9,
    'f1': 0.6206896551724138,
    'accuracy': None}