"""
Classifies: CHEBI:33859 aromatic carboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_aromatic_carboxylic_acid(smiles: str):
    """
    Determines if a molecule is an aromatic carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic carboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"
    
    # Generate the aromatic ring information
    rings = mol.GetRingInfo()
    
    # Check if the carboxylic acid group is bonded to an aromatic ring
    aromatic_rings = []
    for ring in rings.AtomRings():
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if all(atom.GetIsAromatic() for atom in atoms):
            aromatic_rings.append(ring)
    
    if not aromatic_rings:
        return False, "No aromatic rings found"
    
    # Check if any carboxylic acid group is directly bonded to an aromatic ring
    for ring in aromatic_rings:
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in ring and neighbor.GetSymbol() == 'C':
                    for n_neighbor in neighbor.GetNeighbors():
                        if n_neighbor.GetSymbol() == 'O' and n_neighbor.GetIdx() not in ring:
                            if n_neighbor.GetNeighbors()[0].GetSymbol() == 'C' and n_neighbor.GetNeighbors()[0].GetIdx() == neighbor.GetIdx():
                                return True, "Aromatic carboxylic acid found"
    
    return False, "Carboxylic acid group not directly bonded to an aromatic ring"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33859',
                          'name': 'aromatic carboxylic acid',
                          'definition': 'Any carboxylic acid in which the '
                                        'carboxy group is directly bonded to '
                                        'an aromatic ring.',
                          'parents': ['CHEBI:33575', 'CHEBI:33659']},
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
    'num_true_positives': 153,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 19,
    'precision': 0.9935064935064936,
    'recall': 0.8895348837209303,
    'f1': 0.9386503067484663,
    'accuracy': None}