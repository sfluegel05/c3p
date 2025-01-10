"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
"""
Classifies: long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.

    A long-chain fatty acyl-CoA(4-) is a fatty acyl-CoA(4-) arising from deprotonation of the phosphate
    and diphosphate OH groups of any long-chain fatty acyl-CoA; major species at pH 7.3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a long-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the CoA moiety
    # CoA moiety includes the ADP and pantetheine units
    coa_smarts = """
    NCC(=O)NCCS
    """
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Define a SMARTS pattern for the thioester linkage
    # This is the linkage between the fatty acyl chain and the CoA via a thioester bond
    thioester_smarts = "C(=O)SC"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage not found"

    # Find the carbonyl carbon attached to sulfur (thioester)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Get the carbon atom index of the carbonyl carbon in the thioester
    carbonyl_carbon_idx = None
    for match in thioester_matches:
        carbonyl_carbon_idx = match[0]  # The first atom in the pattern is the carbonyl carbon
        break  # Assuming only one thioester linkage

    if carbonyl_carbon_idx is None:
        return False, "Could not find carbonyl carbon in thioester linkage"

    # Traverse the acyl chain attached to the carbonyl carbon
    acyl_chain_atoms = set()
    atoms_to_visit = [carbonyl_carbon_idx]
    visited_atoms = set()
    while atoms_to_visit:
        atom_idx = atoms_to_visit.pop()
        if atom_idx in visited_atoms:
            continue
        visited_atoms.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:  # Carbon atoms only
            acyl_chain_atoms.add(atom_idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                bond = mol.GetBondBetweenAtoms(atom_idx, nbr_idx)
                if bond.GetBondType() != Chem.BondType.SINGLE and bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue
                if nbr.GetAtomicNum() == 6 and nbr_idx not in visited_atoms:
                    atoms_to_visit.append(nbr_idx)
    
    # Exclude the carbonyl carbon
    acyl_chain_length = len(acyl_chain_atoms) - 1
    if acyl_chain_length < 12:
        return False, f"Acyl chain length {acyl_chain_length} is too short for long-chain fatty acid"
    
    return True, "Contains CoA moiety with long-chain fatty acyl thioester linkage"

# Metadata (keeping this as per the example)
__metadata__ = {   'chemical_class': {   'id': None,
                              'name': 'long-chain fatty acyl-CoA(4-)',
                              'definition': 'A fatty acyl-CoA(4-) arising from deprotonation of the phosphate '
                                            'and diphosphate OH groups of any long-chain fatty acyl-CoA; '
                                            'major species at pH 7.3.',
                              'parents': []},
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
        'attempt': 1,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None,
        'num_true_positives': None,
        'num_false_positives': None,
        'num_true_negatives': None,
        'num_false_negatives': None,
        'num_negatives': None,
        'precision': None,
        'recall': None,
        'f1': None,
        'accuracy': None}