"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: epoxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid is a heterocyclic fatty acid containing an epoxide ring.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for long hydrocarbon chain (aliphatic chain)
    carbon_chain_length = 0
    for fragment in Chem.GetMolFrags(mol, asMols=True):
        num_carbons = sum(1 for atom in fragment.GetAtoms() if atom.GetAtomicNum() == 6)
        num_heteroatoms = sum(1 for atom in fragment.GetAtoms() if atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 1)
        # Consider as fatty acid chain if it's mainly carbon atoms
        if num_carbons >= 8 and num_heteroatoms <= 2:
            carbon_chain_length = max(carbon_chain_length, num_carbons)
    if carbon_chain_length < 8:
        return False, "No long hydrocarbon chain found (minimum 8 carbons)"

    # Check for epoxide ring (three-membered cyclic ether)
    epoxide_pattern = Chem.MolFromSmarts("[C;R1]1-[O;R1]-[C;R1]1")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide ring found"

    # Ensure the epoxide ring is part of the hydrocarbon chain
    # Check if the epoxide is connected to aliphatic carbons
    epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)
    epoxide_in_chain = False
    for match in epoxide_matches:
        atom_indices = set(match)
        # Check neighboring atoms to see if they are part of the chain
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and not neighbor.IsInRing():
                    epoxide_in_chain = True
                    break
            if epoxide_in_chain:
                break
        if epoxide_in_chain:
            break
    if not epoxide_in_chain:
        return False, "Epoxide ring not part of hydrocarbon chain"

    return True, "Molecule contains a carboxylic acid group and an epoxide ring within a hydrocarbon chain"


__metadata__ = {   'chemical_class': {   'id': None,
                              'name': 'epoxy fatty acid',
                              'definition': 'A heterocyclic fatty acid containing an epoxide ring as part of its structure.',
                              'parents': ['fatty acid', 'heterocyclic compound']},
        'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8,
                      'max_attempts': 5,
                      'max_positive_instances': None,
                      'max_positive_to_test': None,
                      'max_negative_to_test': None,
                      'max_positive_in_prompt': 50,
                      'max_negative_in_prompt': 20,
                      'max_instances_in_prompt': 100,
                      'test_proportion': 0.1},
        'message': None,
        'attempt': 0,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None,
        'num_true_positives': 0,
        'num_false_positives': 0,
        'num_true_negatives': 0,
        'num_false_negatives': 0,
        'num_negatives': None,
        'precision': None,
        'recall': None,
        'f1': None,
        'accuracy': None}