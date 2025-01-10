"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:60940 short-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    A short-chain fatty acyl-CoA is a fatty acyl-CoA that results from the formal condensation
    of the thiol group of coenzyme A with the carboxy group of any short-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for Coenzyme A moiety (simplified pattern matching key features)
    coA_smarts_str = "NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]1O"
    coA_smarts = Chem.MolFromSmarts(coA_smarts_str)
    if coA_smarts is None:
        return False, "Invalid CoA SMARTS pattern"

    # Check for Coenzyme A substructure
    if not mol.HasSubstructMatch(coA_smarts):
        return False, "Coenzyme A moiety not found"

    # Define SMARTS pattern for thioester linkage (C(=O)S)
    thioester_smarts_str = "C(=O)S"
    thioester_smarts = Chem.MolFromSmarts(thioester_smarts_str)
    if thioester_smarts is None:
        return False, "Invalid thioester SMARTS pattern"

    # Find the thioester linkage
    thioester_matches = mol.GetSubstructMatches(thioester_smarts)
    if not thioester_matches:
        return False, "Thioester linkage not found"

    # For each thioester linkage, attempt to identify the acyl chain
    for match in thioester_matches:
        carbonyl_c_idx = match[0]  # Carbon of C=O
        sulfur_idx = match[2] if len(match) > 2 else match[1]  # S atom index

        # Get the bond between carbonyl carbon and sulfur
        bond = mol.GetBondBetweenAtoms(carbonyl_c_idx, sulfur_idx)
        if bond is None:
            continue

        # Break the bond to separate acyl chain from CoA
        mol_frag = Chem.FragmentOnBonds(mol, [bond.GetIdx()])
        frags = Chem.GetMolFrags(mol_frag, asMols=True, sanitizeFrags=True)

        # Identify the fragment that contains the acyl chain
        acyl_chain = None
        for frag in frags:
            if frag.HasSubstructMatch(thioester_smarts):  # Skip fragment containing thioester
                continue
            # Assume the fragment without phosphorus atoms is the acyl chain
            if not any(atom.GetAtomicNum() == 15 for atom in frag.GetAtoms()):
                acyl_chain = frag
                break

        if acyl_chain is None:
            continue

        # Count the number of carbons in the acyl chain (excluding the carbonyl carbon)
        num_carbons = sum(1 for atom in acyl_chain.GetAtoms() if atom.GetAtomicNum() == 6)

        if 1 <= num_carbons <= 5:
            return True, f"Found short-chain fatty acyl-CoA with acyl chain length {num_carbons}"
        else:
            return False, f"Acyl chain length is {num_carbons}, not short-chain"

    return False, "Acyl chain not identified or not short-chain"

# Add metadata if required
__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:60940',
        'name': 'short-chain fatty acyl-CoA',
        'definition': 'A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any short-chain fatty acid.',
        'parents': ['CHEBI:37554', 'CHEBI:57288']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}