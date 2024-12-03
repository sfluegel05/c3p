"""
Classifies: CHEBI:84135 L-serine derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_L_serine_derivative(smiles: str):
    """
    Determines if a molecule is an L-serine derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an L-serine derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the L-serine substructure
    l_serine_smiles = "N[C@@H](CO)C(=O)O"
    l_serine_mol = Chem.MolFromSmiles(l_serine_smiles)
    
    # Check if the molecule contains the L-serine substructure
    if not mol.HasSubstructMatch(l_serine_mol):
        return False, "Does not contain L-serine substructure"
    
    # Check for modifications at the amino group or carboxy group
    l_serine_matches = mol.GetSubstructMatches(l_serine_mol)
    for match in l_serine_matches:
        amino_group = mol.GetAtomWithIdx(match[0])
        carboxy_group = mol.GetAtomWithIdx(match[2])
        
        # Check for modifications at the amino group
        amino_neighbors = [n.GetSymbol() for n in amino_group.GetNeighbors() if n.GetIdx() not in match]
        if any(n != 'H' for n in amino_neighbors):
            return True, "Modified at the amino group"
        
        # Check for modifications at the carboxy group
        carboxy_neighbors = [n.GetSymbol() for n in carboxy_group.GetNeighbors() if n.GetIdx() not in match]
        if any(n != 'H' for n in carboxy_neighbors):
            return True, "Modified at the carboxy group"
    
    # Check for replacement of any hydrogen by a heteroatom
    heteroatoms = ['N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I']
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'H':
            continue
        neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
        if any(n in heteroatoms for n in neighbors):
            return True, "Hydrogen replaced by a heteroatom"
    
    return False, "No modifications detected"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:84135',
                          'name': 'L-serine derivative',
                          'definition': 'A proteinogenic amino acid derivative '
                                        'resulting from reaction of L-serine '
                                        'at the amino group or the carboxy '
                                        'group, or from the replacement of any '
                                        'hydrogen of L-serine by a heteroatom.',
                          'parents': ['CHEBI:26649', 'CHEBI:83811']},
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
    'num_true_positives': 47,
    'num_false_positives': 9,
    'num_true_negatives': 11,
    'num_false_negatives': 2,
    'precision': 0.8392857142857143,
    'recall': 0.9591836734693877,
    'f1': 0.8952380952380952,
    'accuracy': None}