"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:57618 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.
    A long-chain fatty acyl-CoA(4-) is a fatty acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate OH groups of any long-chain fatty acyl-CoA; major species at pH 7.3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for thioester linkage (S-C(=O)-C)
    thioester_pattern = Chem.MolFromSmarts("S[C](=O)C")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found between fatty acid and CoA"

    # Assume the first thioester match is the acyl linkage
    sulfur_idx = thioester_matches[0][0]
    carbonyl_c_idx = thioester_matches[0][1]
    acyl_c_idx = thioester_matches[0][2]

    # Break the bond between sulfur and carbonyl carbon to separate acyl chain
    bond_to_break = mol.GetBondBetweenAtoms(sulfur_idx, carbonyl_c_idx).GetIdx()
    fragmented_mol = Chem.FragmentOnBonds(mol, [bond_to_break])
    frags = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)

    # Identify acyl chain fragment (without phosphorus atoms)
    acyl_chain_frag = None
    coa_frag = None
    for frag in frags:
        atom_nums = [atom.GetAtomicNum() for atom in frag.GetAtoms()]
        if 15 in atom_nums:
            coa_frag = frag
        else:
            acyl_chain_frag = frag

    if acyl_chain_frag is None or coa_frag is None:
        return False, "Could not separate acyl chain and CoA fragments"

    # Count the number of carbons in acyl chain
    num_carbons = sum(1 for atom in acyl_chain_frag.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 12:
        return False, f"Acyl chain has {num_carbons} carbons, not a long-chain fatty acid"

    # Check for CoA moiety in the molecule
    coa_pattern = Chem.MolFromSmarts("NC(=O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](CO[P](=O)(O)OP(=O)(O)OP(=O)(O)OC[C@H]2O[C@H](n3cnc4c(N)ncnc43)[C@H](O)[C@H]2O)[C@H](O)[C@@H]1O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Check that phosphate and diphosphate groups are deprotonated ([O-])
    deprotonated_phosphate_pattern = Chem.MolFromSmarts("P(=O)([O-])([O-])[O-]")
    phosphate_matches = mol.GetSubstructMatches(deprotonated_phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, "Phosphate groups are not fully deprotonated"

    # Check the total formal charge
    total_charge = Chem.GetFormalCharge(mol)
    if total_charge != -4:
        return False, f"Molecule has net charge {total_charge}, expected -4"

    return True, "Molecule is a long-chain fatty acyl-CoA(4-)"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:57618',
                              'name': 'long-chain fatty acyl-CoA(4-)',
                              'definition': 'A fatty acyl-CoA(4-) arising from the deprotonation of the phosphate and diphosphate OH groups of any long-chain fatty acyl-CoA; major species at pH 7.3.',
                              'parents': ['CHEBI:59474']},
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