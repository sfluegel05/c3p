"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
"""
Classifies: CHEBI:73048 O-acyl-L-carnitine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    An O-acyl-L-carnitine is an O-acylcarnitine with the carnitine component in the L-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the carnitine backbone with L-configuration
    # The pattern should match: [N+](C)(C)C[C@H](CC([O-])=O)OC(=O)
    carnitine_pattern = Chem.MolFromSmarts("[N+X4](C)(C)C[C@H](CC([O-])=O)OC(=O)")
    if not mol.HasSubstructMatch(carnitine_pattern):
        return False, "No L-carnitine backbone found"

    # Check for the presence of a quaternary ammonium group (N+)
    ammonium_pattern = Chem.MolFromSmarts("[N+X4]")
    if not mol.HasSubstructMatch(ammonium_pattern):
        return False, "No quaternary ammonium group found"

    # Check for the carboxylate group (COO-)
    carboxylate_pattern = Chem.MolFromSmarts("CC([O-])=O")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"

    # Check for the ester bond (O-acyl group)
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester bond found"

    # Check for the L-configuration by ensuring the chiral center is correctly configured
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not chiral_centers:
        return False, "No chiral center found"
    
    # Ensure the chiral center is in the L-configuration
    # The chiral center should be the carbon attached to the nitrogen and the carboxylate group
    for center in chiral_centers:
        if center[1] != 'R':  # L-configuration corresponds to 'R' in RDKit
            return False, "Chiral center not in L-configuration"

    return True, "Contains L-carnitine backbone with O-acyl group, quaternary ammonium, and carboxylate group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:73048',
                          'name': 'O-acyl-L-carnitine',
                          'definition': 'An O-acylcarnitine in which the carnitine component has L-configuration.',
                          'parents': ['CHEBI:73047', 'CHEBI:73046']},
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