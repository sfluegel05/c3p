"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: CHEBI:35748 ether lipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid is a lipid with a glycerol backbone where at least one carbon atom is bonded to an alkyl chain via an ether linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for ether linkages (C-O-C)
    ether_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if len(ether_matches) == 0:
        return False, "No ether linkages found"

    # Check for at least one ether linkage connected to the glycerol backbone
    glycerol_atoms = mol.GetSubstructMatch(glycerol_pattern)
    ether_connected_to_glycerol = False
    for match in ether_matches:
        if any(atom in glycerol_atoms for atom in match):
            ether_connected_to_glycerol = True
            break
    if not ether_connected_to_glycerol:
        return False, "Ether linkage not connected to glycerol backbone"

    # Check for long alkyl chains (at least 4 carbons)
    alkyl_chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]")
    alkyl_chain_matches = mol.GetSubstructMatches(alkyl_chain_pattern)
    if len(alkyl_chain_matches) < 1:
        return False, "No long alkyl chains found"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chains too short to be alkyl chains"

    # Check molecular weight - ether lipids typically >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for ether lipid"

    return True, "Contains glycerol backbone with at least one ether linkage to an alkyl chain"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35748',
                          'name': 'ether lipid',
                          'definition': 'A lipid similar in structure to a '
                                        'glycerolipid but in which one or more '
                                        'of the carbon atoms on glycerol is '
                                        'bonded to an alkyl chain via an ether '
                                        'linkage, as opposed to the usual ester '
                                        'linkage.',
                          'parents': ['CHEBI:35747', 'CHEBI:76579']},
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