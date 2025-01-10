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
    # CoA has an ADP (adenosine diphosphate) moiety and a pantetheine moiety connected via phosphate groups
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H](n2cnc3c(N)ncnc23)[C@H](O)[C@@H]1OP(=O)([O-])[O-]")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Define a SMARTS pattern for the thioester linkage
    # This is the linkage between the fatty acyl chain and the CoA via a thioester bond
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage not found"

    # Find the fatty acyl chain attached via thioester linkage
    # Identify the carbonyl carbon attached to sulfur (thioester)
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # Get the carbon atom of the acyl chain
    carbonyl_carbons = [match[0] for match in thioester_matches]
    acyl_lengths = []
    for carbon_idx in carbonyl_carbons:
        # Obtain the acyl chain attached to the carbonyl carbon
        carbon = mol.GetAtomWithIdx(carbon_idx)
        # Traverse the chain away from the carbonyl carbon
        chain_atoms = set()
        atoms_to_visit = [nbr.GetIdx() for nbr in carbon.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != carbon_idx]
        while atoms_to_visit:
            atom_idx = atoms_to_visit.pop()
            if atom_idx not in chain_atoms:
                chain_atoms.add(atom_idx)
                atom = mol.GetAtomWithIdx(atom_idx)
                for nbr in atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr.GetAtomicNum() == 6 and nbr_idx != atom_idx and nbr_idx not in chain_atoms:
                        atoms_to_visit.append(nbr_idx)
        acyl_length = len(chain_atoms)
        acyl_lengths.append(acyl_length)

    if not acyl_lengths:
        return False, "Acyl chain not found"

    # Check if any of the acyl chains are long-chain (12 or more carbons)
    if all(length < 12 for length in acyl_lengths):
        return False, f"Acyl chain length(s) {acyl_lengths} are too short for long-chain fatty acid"

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
        'attempt': 0,
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