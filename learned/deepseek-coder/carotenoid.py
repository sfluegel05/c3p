"""
Classifies: CHEBI:23044 carotenoid
"""
"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    A carotenoid is a tetraterpenoid (C40) derived from psi,psi-carotene, with a long polyene chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a long polyene chain (at least 8 conjugated double bonds)
    polyene_pattern = Chem.MolFromSmarts("[CH]=[CH]-[CH]=[CH]-[CH]=[CH]-[CH]=[CH]-[CH]=[CH]-[CH]=[CH]")
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "No long polyene chain found (at least 8 conjugated double bonds required)"

    # Check molecular weight - carotenoids typically have a molecular weight around 536-600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 650:
        return False, f"Molecular weight {mol_wt:.2f} Da is outside the typical range for carotenoids (500-650 Da)"

    # Count carbons - carotenoids typically have 40 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 35 or c_count > 45:
        return False, f"Found {c_count} carbons, expected around 40 for carotenoids"

    # Check for the presence of cyclic structures (optional, as some carotenoids are acyclic)
    cyclic_pattern = Chem.MolFromSmarts("[R]")
    if not mol.HasSubstructMatch(cyclic_pattern):
        return False, "No cyclic structures found (some carotenoids are cyclic)"

    # Check for functional groups like hydroxyl, carbonyl, etc. (optional)
    functional_group_pattern = Chem.MolFromSmarts("[OH,OX1,C=O]")
    if not mol.HasSubstructMatch(functional_group_pattern):
        return False, "No functional groups found (some carotenoids have hydroxyl or carbonyl groups)"

    return True, "Contains a long polyene chain with conjugated double bonds, typical of carotenoids"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23044',
                          'name': 'carotenoid',
                          'definition': 'One of a class of tetraterpenoids (C40), '
                                        'formally derived from the acyclic parent, '
                                        'psi,psi-carotene by hydrogenation, '
                                        'dehydrogenation, cyclization, oxidation, '
                                        'or combination of these processes. This '
                                        'class includes carotenes, xanthophylls '
                                        'and certain compounds that arise from '
                                        'rearrangement of the skeleton of '
                                        'psi,psi-carotene or by loss of part of '
                                        'this structure. Retinoids are excluded.'},
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
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}